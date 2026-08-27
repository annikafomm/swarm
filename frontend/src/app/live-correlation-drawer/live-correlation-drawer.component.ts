import {
  Component,
  Input,
  Output,
  EventEmitter,
  OnChanges,
  SimpleChanges,
  ViewChild,
  ElementRef,
  OnDestroy
} from '@angular/core';
import { CommonModule } from '@angular/common';
import { HttpClient } from '@angular/common/http';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import * as Plotly from 'plotly.js-dist-min';

export interface SpatialFeatureRef {
  id: string;
  name: string;
}

@Component({
  selector: 'app-live-correlation-drawer',
  standalone: true,
  imports: [
    CommonModule,
    MatIconModule,
    MatButtonModule,
    MatTooltipModule,
    MatProgressSpinnerModule
  ],
  templateUrl: './live-correlation-drawer.component.html',
  styleUrls: ['./live-correlation-drawer.component.scss']
})
export class LiveCorrelationDrawerComponent implements OnChanges, OnDestroy {
  @Input() datasetId: string | undefined;
  @Input() featureA: SpatialFeatureRef | null = null;
  @Input() featureB: SpatialFeatureRef | null = null;
  @Input() isVisible: boolean = false;
  @Input() isBivariateActive: boolean = false;

  @Output() cellSelected = new EventEmitter<string>();
  @Output() applyBivariateColorRequested = new EventEmitter<{ featA: SpatialFeatureRef; featB: SpatialFeatureRef; enabled?: boolean }>();
  @Output() drawerClosed = new EventEmitter<void>();

  @ViewChild('drawerScatterContainer') scatterContainer!: ElementRef<HTMLDivElement>;

  public isCollapsed: boolean = false;
  public isLoading: boolean = false;
  public scatterData: any = null;
  public errorMsg: string = '';

  private currentFetchSub: any = null;

  constructor(private http: HttpClient) {}

  public onClose(): void {
    this.drawerClosed.emit();
  }

  public toggleBivariateMap(): void {
    if (!this.featureA || !this.featureB) return;
    this.isBivariateActive = !this.isBivariateActive;
    this.applyBivariateColorRequested.emit({
      featA: this.featureA,
      featB: this.featureB,
      enabled: this.isBivariateActive
    });
  }

  ngOnChanges(changes: SimpleChanges): void {
    const prevA = changes['featureA']?.previousValue?.id;
    const currA = changes['featureA']?.currentValue?.id;
    const prevB = changes['featureB']?.previousValue?.id;
    const currB = changes['featureB']?.currentValue?.id;
    const prevDs = changes['datasetId']?.previousValue;
    const currDs = changes['datasetId']?.currentValue;

    const featAChanged = prevA !== currA && !!currA;
    const featBChanged = prevB !== currB && !!currB;
    const datasetChanged = prevDs !== currDs && !!currDs;
    const visChanged = changes['isVisible'] && changes['isVisible'].currentValue && !changes['isVisible'].previousValue;

    if (this.isVisible && (featAChanged || featBChanged || datasetChanged || visChanged)) {
      this.loadCorrelationData();
    }
  }

  ngOnDestroy(): void {
    if (this.scatterContainer?.nativeElement) {
      Plotly.purge(this.scatterContainer.nativeElement);
    }
  }

  public toggleCollapse(): void {
    this.isCollapsed = !this.isCollapsed;
    if (!this.isCollapsed && this.scatterData) {
      setTimeout(() => this.renderScatterplot(), 50);
    }
  }

  public loadCorrelationData(): void {
    if (!this.datasetId || !this.featureA?.id || !this.featureB?.id) {
      return;
    }

    this.isLoading = true;
    this.errorMsg = '';

    const url = `/api/datasets/${encodeURIComponent(this.datasetId)}/spatial_correlation_pair?feature_id_a=${encodeURIComponent(this.featureA.id)}&feature_id_b=${encodeURIComponent(this.featureB.id)}`;

    if (this.currentFetchSub) {
      this.currentFetchSub.unsubscribe();
    }

    this.currentFetchSub = this.http.get<any>(url, { withCredentials: true }).subscribe({
      next: (res) => {
        this.isLoading = false;
        this.scatterData = res;
        if (!this.isCollapsed) {
          setTimeout(() => this.renderScatterplot(), 50);
        }
      },
      error: (err) => {
        this.isLoading = false;
        this.errorMsg = 'Failed to calculate cross-map correlation';
        console.error('Failed to load drawer scatter data:', err);
      }
    });
  }

  private renderScatterplot(): void {
    if (!this.scatterContainer?.nativeElement || !this.scatterData) return;

    const container = this.scatterContainer.nativeElement;
    const points = this.scatterData.points || [];
    const stats = this.scatterData.stats || {};

    const xVals = points.map((p: any) => p.x);
    const yVals = points.map((p: any) => p.y);
    const cellIds = points.map((p: any) => p.cell_id);
    const clusters = points.map((p: any) => p.cluster);

    const hoverTexts = points.map(
      (p: any) =>
        `<b>Cell:</b> ${p.cell_id}<br>` +
        `<b>${this.featureA?.name}:</b> ${typeof p.x === 'number' ? p.x.toFixed(3) : p.x}<br>` +
        `<b>${this.featureB?.name}:</b> ${typeof p.y === 'number' ? p.y.toFixed(3) : p.y}`
    );

    const data: Partial<Plotly.PlotData>[] = [
      {
        type: 'scatter',
        mode: 'markers',
        x: xVals,
        y: yVals,
        text: hoverTexts,
        customdata: cellIds,
        hoverinfo: 'text',
        marker: {
          size: 5,
          color: clusters as any,
          colorscale: 'Viridis',
          opacity: 0.75
        },
        name: 'Cells'
      }
    ];

    if (points.length > 1 && stats.slope !== undefined && stats.intercept !== undefined) {
      const minX = Math.min(...xVals);
      const maxX = Math.max(...xVals);
      data.push({
        type: 'scatter',
        mode: 'lines',
        x: [minX, maxX],
        y: [minX * stats.slope + stats.intercept, maxX * stats.slope + stats.intercept],
        line: { color: '#ef4444', width: 2, dash: 'dot' },
        name: 'OLS Trend'
      });
    }

    const layout: Partial<Plotly.Layout> = {
      margin: { t: 15, l: 45, r: 15, b: 35 },
      height: 220,
      showlegend: false,
      xaxis: {
        title: { text: this.featureA?.name || 'Map 1 (Left)', font: { size: 10 } },
        automargin: true,
        tickfont: { size: 9 }
      },
      yaxis: {
        title: { text: this.featureB?.name || 'Map 2 (Right)', font: { size: 10 } },
        automargin: true,
        tickfont: { size: 9 }
      }
    };

    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, {
      responsive: true,
      displayModeBar: false
    }).then(() => {
      (container as any).on('plotly_click', (eventData: any) => {
        if (eventData && eventData.points && eventData.points.length > 0) {
          const pt = eventData.points[0];
          if (pt.customdata) {
            this.cellSelected.emit(String(pt.customdata));
          }
        }
      });
    });
  }

  public onApplyBivariate(): void {
    this.toggleBivariateMap();
  }
}

