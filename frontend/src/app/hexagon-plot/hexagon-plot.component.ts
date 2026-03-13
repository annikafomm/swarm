import { Component, OnInit, OnDestroy, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { KeyValue } from '@angular/common';
import { Data, Router } from '@angular/router';
import { map, Observable, Subject, Subscription, takeUntil } from 'rxjs';

// Visualizations
import * as d3 from 'd3';
import * as Plotly from 'plotly.js-dist-min';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import Shepherd from 'shepherd.js';

// Services
import { DatasetService } from '../datasets.service';
import { Dataset } from '../datasets.service';
import { SessionService } from '../session.service';
import { GeoDataService } from '../geo-data.service';
import { TranslationService } from '../translation.service';
import { TranslatePipe } from '../translate.pipe';
import { PathsService } from '../paths.service';
import { DEFAULT_PATHS } from '../constants';

// Material
import { MatButtonModule } from '@angular/material/button';
import { MatIconModule } from '@angular/material/icon';
import { MatTooltipModule } from '@angular/material/tooltip';
import { MatDialogModule } from '@angular/material/dialog';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatFormField } from '@angular/material/form-field';
import { MatOptgroup, MatOption } from '@angular/material/autocomplete';
import { MatLabel } from '@angular/material/form-field';
import { MatSelect } from '@angular/material/select';
import { MatExpansionModule } from '@angular/material/expansion';
import { MatTableModule } from '@angular/material/table';
import { MatDividerModule } from '@angular/material/divider';
import { MatTabsModule, MatTabHeader } from '@angular/material/tabs';
import { MatTabChangeEvent, MatTabGroup } from '@angular/material/tabs';
import { DomSanitizer, SafeResourceUrl } from '@angular/platform-browser';
import { MatInputModule } from '@angular/material/input';



@Component({
  selector: 'app-hexagon-plot',
  imports: [CommonModule, FormsModule, FilterableTableComponent, TranslatePipe, MatButtonModule, MatIconModule, MatTooltipModule, MatDialogModule, MatProgressSpinnerModule, MatOptgroup, MatFormField, MatLabel, MatOption, MatSelect, MatExpansionModule, MatTableModule, MatDividerModule, MatTabsModule, MatTabHeader, MatInputModule],
  standalone: true,
  templateUrl: './hexagon-plot.component.html',
  styleUrls: ['./hexagon-plot.component.scss'],
})
export class HexagonPlotComponent implements OnInit, OnDestroy, AfterViewInit {
  @ViewChild('aucell_graph_genie3', { static: false }) aucellGraphGenie3Element!: ElementRef<HTMLElement>;
  @ViewChild('aucell_graph_sponge', { static: false }) aucellGraphSpongeElement!: ElementRef<HTMLElement>;
  @ViewChild(MatTabGroup, { static: false }) tabGroup?: MatTabGroup;
  @ViewChild('dgeaHeatmap', { static: false }) dgeaHeatmapElement!: ElementRef<HTMLElement>;
  private _resizeHandler: any = null;
  private sub!: Subscription;
  footprintPlotUrls: SafeResourceUrl[] = [];
  onDemandFootprintUrls: SafeResourceUrl[] = [];
  availableMotifs: string[] = [];
  footprintMotif: string = '';
  footprintClusterBy: string = 'cell_type';
  isComputingFootprint: boolean = false;
  footprintComputeError: string = '';

  builtinDatasets$: Observable<Dataset[]>;
  uploadedDatasets$: Observable<Dataset[]>;
  tangram_Datasets$: Observable<Dataset[]>;
  selectedDataset: Dataset | null = null;
  selectedDatasetCompare: Dataset | null = null;

  compareMode: boolean = false;;

  // Track component destruction to avoid setting state on an unmounted component
  private destroy$ = new Subject<void>();

  constructor(
    private router: Router,
    private http: HttpClient,
    private sessionService: SessionService,
    private geoDataService: GeoDataService,
    private translationService: TranslationService,
    private datasetService: DatasetService,
    private pathsService: PathsService,
    private sanitizer: DomSanitizer,
  ) {

    // Setup dataset observables
    this.builtinDatasets$ = this.datasetService.availableDatasets$.pipe(
      map(datasets => datasets.filter(d => d.type === 'builtin'))
    );
    this.uploadedDatasets$ = this.datasetService.availableDatasets$.pipe(
      map(datasets => datasets.filter(d => d.type === 'uploaded'))
    );
    this.tangram_Datasets$ = this.datasetService.availableDatasets$.pipe(
      map(datasets => datasets.filter(d => d.tangram_adata_path))
    );

  }


  // Define to use Math functions in the html template
  public Math = Math;

  // GeoJson
  public dataPath = DEFAULT_PATHS.hexagonPath;
  public dataSetTitle =
    this.dataPath.split('/').pop()?.replace('.geojson', '') || 'Hexagon Plot';

  // Map svg and g elements
  private svg!: d3.Selection<SVGSVGElement, any, any, any>;
  private g!: d3.Selection<SVGGElement, any, any, any>;
  private g_compare!: d3.Selection<SVGGElement, any, any, any>;
  private svg_compare!: d3.Selection<SVGSVGElement, any, any, any>;

    // Nested g elements that contain the actual paths
    private g_paths!: d3.Selection<SVGGElement, any, any, any>;
    private g_paths_compare!: d3.Selection<SVGGElement, any, any, any>;
    // ======= Xenium performance state =======
    private fullFeatures: CellFeature[] = [];
    private isXenium = false;

  private baseLayer!: d3.Selection<SVGGElement, null, any, any>;
  private detailLayer!: d3.Selection<SVGGElement, null, any, any>;

  private currentTransform = d3.zoomIdentity;
  private currentPathGenerator!: d3.GeoPath<any, CellFeature>;

  private detailSize = 80;

  private detailVisible = false;
  private detailScreenPos: { x: number; y: number } | null = null;
  private keydownHandler?: (event: KeyboardEvent) => void;
  // =======================================


    public selectedCell: CellFeature | null = null;
    public selectedCellCompare: CellFeature | null = null;
    public selectedCluster: number | null = null;
    public colorByProperty = 'regulatory_scores';
    public selectedGeneSetGenie3: string | null = null;
    public selectedGeneSetSponge: string | null = null;
    public selectedRegulatoryScore: string | null = null;

  // Data sources for the two tables
  public genie3RawData: TableData = {};
  public spongeRawData: TableData = {};

  // Column lists for the two tables
  public genie3Elements: string[] = [];
  public spongeElements: string[] = [];


  private previousGeneSetGenie3: string | null = null;
  private previousGeneSetSponge: string | null = null;
  private requestTokens: { [key: string]: number } = {};
  public dgeaReady: boolean = false;

  public selectedInterval: number = 0;
  public features: CellFeature[] = []; // public so that filterable table can update it
  public meta: { [key: string]: any } = {};

    // Selected groups for the DGEA comparison (bound to the dropdowns)
    public selectedDgeaGroup1: string | null = null;
    public selectedDgeaGroup2: string | null = null;

    private hiddenPropKeys = new Set<string>([]);

  // Selected groups for the DGEA comparison (bound to the dropdowns)
  getSelectedDgeaHeatmap(): any | null {
    const cmp = this.getSelectedDgeaComparison();
    if (!cmp || cmp['skipped']) return null;
    return cmp['heatmap_context'] ?? null;
  }

  // Returns all available cell type levels for the dropdown selectors
  getDgeaCellTypeLevels(): string[] {
    return this.meta?.['dgea']?.['cell_type']?.['levels'] ?? [];
  }

  // Returns the map of all DGEA comparisons
  getDgeaComparisonMap(): { [key: string]: any } {
    return this.meta?.['dgea']?.['cell_type']?.['comparisons'] ?? {};
  }

    // Re-render the DGEA heatmap when the user changes the group selections
    public onDgeaSelectionChange(): void {
        setTimeout(() => this.renderDgeaHeatmap(), 0);
    }

    private getLeidenClusterAnnotation(clusterId: number | null | undefined): any | null {
      if (clusterId === null || clusterId === undefined) return null;
      return this.meta?.['leiden_cluster_annotations']?.[String(clusterId)] ?? null;
    }

  // Create comparison ID matching the backend JSON format
  private makeComparisonId(group1: string, group2: string): string {
    const safe = (x: string) => x.replace(/[^A-Za-z0-9]+/g, '_');
    return `${safe(group1)}__vs__${safe(group2)}`;
  }

  // Render the context heatmap
  private renderDgeaHeatmap(): void {
    const container = this.dgeaHeatmapElement?.nativeElement;
    if (!container) return;

    const hm = this.getSelectedDgeaHeatmap();

    if (!hm || !hm.groups || !hm.rows || !hm.rows.length) {
      Plotly.purge(container);
      return;
    }

    const x = hm.groups as string[];
    const y = hm.rows.map((r: any) => r.gene);
    const z = hm.rows.map((r: any) => r.scaled);
    const raw = hm.rows.map((r: any) => r.raw);

    const data: Partial<Plotly.PlotData>[] = [
      {
        type: 'heatmap',
        x,
        y,
        z,
        customdata: raw,
        colorscale: 'RdBu',
        reversescale: true,
        hovertemplate:
          'Gene: %{y}<br>' +
          'Group: %{x}<br>' +
          'Scaled expression: %{z:.2f}<br>' +
          'Mean expression: %{customdata:.2f}<extra></extra>',
      }
    ];

    const layout: Partial<Plotly.Layout> = {
      margin: { t: 30, l: 140, r: 20, b: 100 },
      height: Math.max(420, y.length * 22),
      xaxis: {
        title: { text: 'Cell type' },
        tickangle: -45,
        automargin: true
      },
      yaxis: {
        title: { text: 'Genes' },
        automargin: true,
        autorange: 'reversed'
      }
    };

    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, {
      responsive: true,
      displayModeBar: false
    });
  }

  // Returns all available cell type levels for the dropdown selectors
  getSelectedDgeaComparison(): any | null {
    if (!this.selectedDgeaGroup1 || !this.selectedDgeaGroup2) return null;
    if (this.selectedDgeaGroup1 === this.selectedDgeaGroup2) return null;

    const comps = this.getDgeaComparisonMap();

    const directId = this.makeComparisonId(this.selectedDgeaGroup1, this.selectedDgeaGroup2);
    if (comps[directId]) return comps[directId];

    const reverseId = this.makeComparisonId(this.selectedDgeaGroup2, this.selectedDgeaGroup1);
    if (comps[reverseId]) return comps[reverseId];

    return null;
  }

  // Initialize default selections for the DGEA comparison dropdowns
  initDgeaSelection(): void {
    const levels = this.getDgeaCellTypeLevels();
    if (!levels.length) return;

    if (!this.selectedDgeaGroup1) {
      this.selectedDgeaGroup1 = levels[0];
    }

    if (!this.selectedDgeaGroup2) {
      const firstDifferent = levels.find(x => x !== this.selectedDgeaGroup1);
      this.selectedDgeaGroup2 = firstDifferent ?? null;
    }
  }


  public clusterCells: CellFeature[] = [];
  public clusterCellTypes: {
    type: string;
    count: number;
    percentage: string;
  }[] = [];
  public clusterCentralityAvg: {
    degree_centrality: number;
    average_clustering: number;
    closeness_centrality: number;
  } = {
      degree_centrality: 0,
      average_clustering: 0,
      closeness_centrality: 0,
    };

  // GRN network and genesets

  public genie3Network: genie3RegGraphConnection[] = [];
  public spongeNetwork: spongeRegGraphConnection[] = [];
  public geneSetsGenie3: { [regulator: string]: string[] } = {};
  public geneSetsSponge: { [regulator: string]: string[] } = {};
  public genie3Width: number = 600;
  public spongeWidth: number = 600;

  // Slider params
  public genie3WeightCutoff: number = 0.5;
  public spongePValueCutoff: number = 0.05;
  public genie3MinEdges: number = 25;
  public spongeMinEdges: number = 25;
  public genie3SliderData: { step: number; min_border: number; max_border: number; default_value: number } | null = null;
  public spongeSliderData: { step: number; min_border: number; max_border: number; default_value: number } | null = null;

  // Loading screen trackers
  public isLoadingHexagons: boolean = true;
  public isLoadingCompare: boolean = false;
  public isLoadingSponge: boolean = false;
  public isLoadingGenie3: boolean = false;

  // Co-occurrence table
  public coOccurrenceData: number[] = [];
  public coOccurrenceColumns: string[] = [];
  public coOccurrenceThreshold: number = 0.5;
  public maxInterval: number = 49;
  public clusterCount: number = 10;

  public colorableProperties: string[] = [
    'cell_type',
    'leiden',
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];

  // keep track of which properties actually have any data available; used to
  // disable/grey‑out menu items and tabs instead of hiding them abruptly.
  public propertyAvailability: { [prop: string]: boolean } = {};

  public leidenCentralityProps = [
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];
  public groupedProperties: { key: string; value: string[] }[] | null = null;
  public ligandReceptorScores: {
    [col: string]: { [index: string]: string | number };
  } | null = null;
  public showGlobalLigandReceptorScores: boolean = true;
  public showMoranI: boolean = true;
  public showChromvarMoranI: boolean = true;
  public colorScale = d3
    .scaleOrdinal<string>()
    .range([
      '#ff9800',
      '#1976d2',
      '#ff6f00',
      '#00bcd4',
      '#ffa726',
      '#9c27b0',
      '#ff8a65',
      '#4caf50',
    ]);
  private continuousColorScale = d3.scaleSequential(d3.interpolateYlOrRd);
  // Separate scales for the compare view to avoid cross-contamination
  public colorScaleCompare = d3.scaleOrdinal<string>().range([
    '#ff9800',
    '#1976d2',
    '#ff6f00',
    '#00bcd4',
    '#ffa726',
    '#9c27b0',
    '#ff8a65',
    '#4caf50',
  ]);
  // Yellow continuous color palette
  private continuousColorScaleCompare = d3.scaleSequential(d3.interpolateYlOrRd);
  public currentLegendDomain: any[] = [];
  public currentLegendType: 'continuous' | 'categorical' = 'categorical';

  // Comparison
  public selectedCompareView: string = 'regulatory_scores';
  public currentCompareLegendType: 'continuous' | 'categorical' = 'categorical';
  public currentLegendDomainCompare: any[] = [];
  private dataCompare: GeoJsonData | null = null;

  // allow nested tables like for differential motif activity view
  public asTableData(value: unknown): {
    [col: string]: { [index: string]: string | number }
  } | string[] {
    return value as { [col: string]: { [index: string]: string | number } } | string[];
  }

  ngOnInit(): void {
    // Initialize with default builtin dataset if no dataset is selected
    this.datasetService.availableDatasets$
      .pipe(
        takeUntil(this.destroy$)
      )
      .subscribe(datasets => {
        console.log('Available datasets:', datasets);
        // Auto-select the first builtin dataset if none selected
        if (datasets.length > 0 && !this.selectedDataset) {
          const builtinDataset = datasets.find(d => d.type === 'builtin');
          if (builtinDataset) {
            this.selectedDataset = builtinDataset;
            this.datasetService.selectDataset(builtinDataset);
          }
        }
      });

    this.datasetService.selectedDataset$
      .pipe(takeUntil(this.destroy$))
      .subscribe(dataset => {
        this.selectedDataset = dataset;
        if (dataset) {
          this.updatePathsFromDataset(dataset);
        }
      });

    this.datasetService.selectedDatasetCompare$
      .pipe(takeUntil(this.destroy$))
      .subscribe(dataset => {
        this.selectedDatasetCompare = dataset;
        if (dataset && this.compareMode) {
          this.reloadComparisonView();
        }
      });

    this.isLoadingHexagons = true;
    // Subscribe to path changes
    this.sub = this.pathsService.paths$.subscribe(paths => {
      console.log("Paths updated:", paths)
      const hexagonPath = paths.hexagonPath || DEFAULT_PATHS.hexagonPath;

      if (hexagonPath) {
        // Update the component's dataPath
        this.dataPath = hexagonPath;
        console.log('✓ Loading hexagon data from', this.dataPath);

        // Clear hexagon-plot container
        d3.select('#hexbin').selectAll('svg').remove();
        d3.select('#hexbin-compare').selectAll('svg').remove();

        this.footprintPlotUrls = [];
        this.onDemandFootprintUrls = [];
        this.availableMotifs = [];
        this.footprintMotif = '';
        this.footprintComputeError = '';

        // Load and render new data
        this.createHexagonPlot();
        this.loadAndRenderData(this.dataPath);

        this.renderFootprintPlots(this.selectedDataset);
      } else {
        console.warn('✗ No hexagon path available');
      }

    });
  }

  ngOnDestroy(): void {

    this.destroy$.next();
    this.destroy$.complete();

    if (this.sub) {
      this.sub.unsubscribe();
    }
    // remove resize listener if added
    try {
      if (this._resizeHandler) window.removeEventListener('resize', this._resizeHandler as any);
    } catch (e) { }

    try {
      if (this.keydownHandler) {
        window.removeEventListener('keydown', this.keydownHandler);
      }
    } catch (e) { }
  }

  ngAfterViewInit(): void {
    this.genie3Width = this.aucellGraphGenie3Element.nativeElement.clientWidth as number;
    this.spongeWidth = this.aucellGraphSpongeElement.nativeElement.clientWidth as number;
    // Force Material Tabs to recalc pagination so arrows appear when needed
    setTimeout(() => this.updateTabPagination(), 50);
    // update on window resize as well
    this._resizeHandler = () => this.updateTabPagination();
    window.addEventListener('resize', this._resizeHandler);
  }

  onDatasetSelected(dataset: Dataset | null): void {
    if (!dataset) return;

    this.selectedDataset = dataset;
    this.datasetService.selectDataset(dataset);  // Networks load automatically
    this.updatePathsFromDataset(dataset);
    this.reloadHexagons();
  }

  // Handle tangram dataset selection - use tangram_adata_path if available
  onDatasetTangramSelected(dataset: Dataset | null): void {
    if (!dataset || !dataset.tangram_adata_path) return;

    // Create a modified dataset with tangram path
    const tangramDataset: Dataset = {
      ...dataset,
      adata_path: dataset.tangram_adata_path
    };

    this.selectedDataset = tangramDataset;
    this.datasetService.selectDataset(tangramDataset);
    this.updatePathsFromDataset(tangramDataset);
    this.reloadHexagons();
  }

  onDatasetCompareSelected(dataset: Dataset | null): void {
    this.selectedDatasetCompare = dataset;
    this.datasetService.selectDatasetCompare(dataset);  // Networks load automatically
    if (dataset) {
      this.reloadComparisonView();
    }
  }

  reloadHexagons(): void {
    if (this.dataPath) {
      this.isLoadingHexagons = true;
      // Clear existing hexagons
      d3.select('#hexbin').selectAll('svg').remove();
      // Load and render new data
      this.createHexagonPlot();
      this.loadAndRenderData(this.dataPath);
    }
  }

  reloadComparisonView(): void {
    if (this.selectedDatasetCompare && this.compareMode) {
      this.isLoadingCompare = true;
      // Clear existing compare hexagons
      d3.select('#hexbin-compare').selectAll('svg').remove();
      // Load and render new compare data
      this.initCompareHexagonPlot();
    }
  }

  private updatePathsFromDataset(dataset: Dataset): void {
    this.pathsService.updatePaths({
      adataMainPath: dataset.adata_path,
      adataComparePath: dataset.adata_path,
      genieFiltPath: dataset.genie_network_path,
      spongeFiltPath: dataset.sponge_network_path,
      hexagonPath: dataset.geojson_path,
    });
  }


  // Track Tab pagination state to change to view accordingly
  private updateTabPagination(): void {
    try {
      if (this.tabGroup && (this.tabGroup as any)._header && typeof (this.tabGroup as any)._header.updatePagination === 'function') {
        (this.tabGroup as any)._header.updatePagination();
      } else if (this.tabGroup && typeof (this.tabGroup as any)._updatePagination === 'function') {
        (this.tabGroup as any)._updatePagination();
      }
    } catch (e) {
      // ignore
    }
  }


  private nextRequestToken(graphType: string): number {
    if (!this.requestTokens[graphType]) this.requestTokens[graphType] = 0;
    return ++this.requestTokens[graphType];
  }

  public openInfoPage(fragmentId?: string): void {
    const frag = fragmentId ? `#${fragmentId}` : '';
    const targetUrl = `/info${frag}`;

    // Navigate via the router (keeps SPA behavior) then force a full-page load
    // to ensure the browser jumps to the fragment reliably.
    this.router.navigate(['/info'], fragmentId ? { fragment: fragmentId } : {}).then(() => {
      // Use replaceState to avoid adding a duplicate history entry, then reload.
      try {
        history.replaceState(null, '', targetUrl);
      } catch { }
      // Force reload so the browser will honor the fragment scroll.
      window.location.href = targetUrl;
    });
  }

  public onTabColorChange(newView: string): void {
    if (this.colorByProperty !== newView) {
      this.colorByProperty = newView;
      this.onColorbyPropertyChange();
    } else {
      // If the colorByProperty is already set, simply force an update.
      this.updateHexColors();
    }
  }

  public onTabChange(event: MatTabChangeEvent): void {
    let newView: string | null = null;
    const tabLabel = event.tab.textLabel;

    // Handle info tabs that need async rendering
    if (tabLabel === 'Cluster Information') {
      // Render neighborhood enrichment when cluster info tab is opened
      // Use longer timeout to ensure Angular has rendered the template
      setTimeout(() => this.renderNhoodHeatmap(), 300);
      return;
    }

    if (tabLabel === 'Cell Information') {
      // Render neighborhood enrichment when cell info tab is opened
      // Use longer timeout to ensure Angular has rendered the template
      setTimeout(() => this.renderNhoodHeatmap(), 300);
      return;
    }

    if (tabLabel === 'DGEA') {
      this.dgeaReady = !!this.meta?.['dgea']?.['cell_type'];

      if (this.dgeaReady) {
        this.initDgeaSelection();
        setTimeout(() => this.renderDgeaHeatmap(), 100);
      }

      return;
    }

    // Handle visualization tabs that change the color property
    switch (tabLabel) {
      case 'Regulatory Scores':
        // Check if regulatory scores data is available
        if (!this.hasRegulatoryScoresData()) {
          console.warn('[Tab Change] Regulatory Scores tab selected but no data available. Skipping view change.');
          return;
        }
        newView = 'regulatory_scores';
        break;
      case 'Co-occurence':
        newView = 'leiden';
        break;
      case 'Gene Expression':
        newView = 'gene_expression';
        break;
      case 'TF Activity':
        newView = 'tf_activity';
        break;
      case 'Pathway Activity':
        newView = 'pathway_activity';
        break;
      case 'Cell Composition TF Activity':
        newView = 'cell_comp_tf_activity_similarity';
        break;
      case 'Ligand-Receptor Relationships':
        newView = 'ligand_receptor_relationships';
        break;
      case "ChromVar spatial correlation : Moran's I / Geary's C":
        newView = 'chromvar_total_sum';
        break;
      case 'Differential Motif Activity':
        // view should be cell type
        newView = 'cell_type';
        break;
      case 'Footprints':
        newView = 'cell_type';
        break;
    }

    console.log('[Tab Change] newView=', newView);

    if (newView && this.colorByProperty !== newView) {
      this.colorByProperty = newView;
      this.onColorbyPropertyChange();
    }
  }

  public onCompareTabChange(event: MatTabChangeEvent): void {
    let newView: string | null = null;

    const label = event.tab.textLabel;

    if (label.startsWith('Compare - ')) {
      const viewName = label.replace('Compare - ', '').toLowerCase().replace(/\s+/g, '_');

      switch (viewName) {
        case 'regulatory_scores':
          // Check if regulatory scores data is available
          if (!this.hasRegulatoryScoresData()) {
            console.warn('[Compare Tab Change] Regulatory Scores tab selected but no data available. Skipping view change.');
            return;
          }
          newView = 'regulatory_scores';
          break;
        case 'gene_expression':
          newView = 'gene_expression';
          break;
        case 'tf_activity':
          newView = 'tf_activity';
          break;
        case 'pathway_activity':
          newView = 'pathway_activity';
          break;
      }
    }

    console.log('[Compare Tab Change] newView=', newView);

    if (newView && this.selectedCompareView !== newView) {
      this.selectedCompareView = newView;
      this.onColorbyPropertyChange();
    }
  }

  /**
   * Check if regulatory scores data is available
   */
  private hasRegulatoryScoresData(): boolean {
    if (!this.meta) return false;

    const hasGenie3 = this.meta['global_regulatory_scores_genie3'] &&
      Object.keys(this.meta['global_regulatory_scores_genie3']).length > 0;
    const hasSponge = this.meta['global_regulatory_scores_sponge'] &&
      Object.keys(this.meta['global_regulatory_scores_sponge']).length > 0;

    return hasGenie3 || hasSponge;
  }

  /**
   * Returns true if any data exists for the given property/column.  This is used
   * by the template to disable tabs and dropdown options when information
   * hasn't been computed or loaded yet.
   */
  public propertyAvailable(prop: string): boolean {
    if (prop === 'co_occurrence') {
      const clusters = this.meta?.['leiden_cluster_annotations'];
      if (!clusters) return false;

      return Object.values(clusters).some(
        (c: any) => Array.isArray(c?.co_occurrence) && c.co_occurrence.length > 0
      );
    }
    return !!this.propertyAvailability[prop];
  }

  /**
   * Compatibility wrapper for the existing regulatory scores tab.  Left in place
   * simply so that the template change is minimal; it now defers to the generic
   * availability check.
   */
  public canViewRegulatoryScores(): boolean {
    return this.propertyAvailable('regulatory_scores');
  }

  private createHexagonPlot(containerName?: string): void {
    if (!containerName) containerName = '#hexbin';

    const width = 500;
    const height = 400;

    // use data-join to reuse an existing svg in the container or create a new one
    const container = d3.select(containerName);
    const svgSel = container
      .selectAll('svg')
      .data([0])
      .join('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', [0, 0, 1200, 1000] as [number, number, number, number])
      .style('background-color', 'white')
      .style('overflow', 'hidden');

    // create (or reuse) a root group inside the svg
    const gSel = svgSel
      .selectAll<SVGGElement, number>('g.root-group')
      .data([0])
      .join('g')
      .attr('class', 'root-group');

    if (containerName === '#hexbin-compare') {
      // Compare view uses its own svg/group
      this.svg_compare = svgSel as unknown as d3.Selection<SVGSVGElement, any, any, any>;
      this.g_compare = gSel as any;
    } else {
      // Main view
      this.svg = svgSel as unknown as d3.Selection<SVGSVGElement, any, any, any>;
      this.g = gSel as any;

      // Xenium layers live only in main view
      this.baseLayer = this.g.selectAll<SVGGElement, unknown>('g.base-layer')
        .data([null])
        .join('g')
        .attr('class', 'base-layer');

      this.detailLayer = this.svg.selectAll<SVGGElement, unknown>('g.detail-layer')
        .data([null])
        .join('g')
        .attr('class', 'detail-layer')
        .attr('clip-path', 'url(#detail-clip)');
    }

    // attach a zoom handler once (reusing svgSel is safe)

    const zoomBehavior = d3
      .zoom<SVGSVGElement, unknown>()
      .scaleExtent([1, 5])
      .extent([
        [0, 0],
        [width, height],
      ])
      .on('zoom', (event) => {
        if (containerName === '#hexbin-compare' && this.g_compare) {
          this.g_compare.attr('transform', event.transform.toString());
        } else if (this.g) {
          this.g.attr('transform', event.transform.toString());
        }
        this.currentTransform = event.transform;
      });

    (svgSel as any).call(zoomBehavior);
  }

  private initCompareHexagonPlot(): void {

    // Do not recurse — if the compare container isn't in the DOM, abort and let the caller retry once.
    const node = document.getElementById('hexbin-compare');
    if (!node) {
      console.warn('hexbin-compare not present in DOM, skipping compare init');
      return;
    }

    // prepare svg/group

    this.createHexagonPlot('#hexbin-compare');

    this.dataCompare = this.geoDataService.getData();
    const featuresToDraw = (this.dataCompare && Array.isArray(this.dataCompare.features)) ? this.dataCompare.features : [];

    if (!this.selectedCompareView) {
      console.log('No compare view selected; skipping compare render');
      return;
    }
    const compareIsContinuous = this.isContinuousScale(this.selectedCompareView, featuresToDraw);
    console.log('[compare init] selectedCompareView=', this.selectedCompareView, 'featuresToDraw=', featuresToDraw.length, 'isContinuous=', compareIsContinuous);
    this.currentCompareLegendType = compareIsContinuous ? 'continuous' : 'categorical';

    const width = 1200;
    const height = 1000;
    const projection = d3.geoIdentity().fitSize([width, height], {
      type: 'FeatureCollection',
      features: featuresToDraw,
    });
    const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

    if (!this.g_compare) {
      console.warn('g_compare not initialized after createHexagonPlot');
      return;
    }

    // draw using the intended feature set (featuresToDraw)
    this.g_paths_compare = this.g_compare.selectAll('g.compare-layer').data([0]).join('g').attr('class', 'compare-layer') as any;

    this.g_paths_compare
      .selectAll('path')
      .data(featuresToDraw)
      .join('path')
      .attr('d', (d) => pathGenerator(d) || '')
      .attr('fill', (d) => {
        const value = (this.leidenCentralityProps.includes(this.selectedCompareView))
          ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[this.selectedCompareView]
          : d.properties[this.selectedCompareView];
        if (this.currentCompareLegendType === 'categorical') {
          return this.colorScaleCompare(String(value));
        } else {
          const num = this.toNumber(value);
          return Number.isFinite(num) ? this.continuousColorScaleCompare(num) : '#ccc';
        }
      })
      .style('opacity', 0.8)
      .on('mouseover', (event, d) => this.mouseOver(event, d))
      .on('mouseleave', (event, d) => this.mouseLeave(event, d))
      .on('click', (event, d) => this.openSidenavCompare(event, d));

    setTimeout(() => { this.isLoadingCompare = false; }, 100);

    setTimeout(() => { this.updateHexColors('#hexbin-compare') }, 0);


  }

  public onCompareMode(): void {

    this.compareMode = !this.compareMode;

    if (this.compareMode) {
      this.isLoadingCompare = true;
      // schedule a single init attempt so Angular has time to render the compare container
      setTimeout(() => this.initCompareHexagonPlot(), 50);


    } else {
      // remove compare svg and clear references
      d3.select('#hexbin-compare').selectAll('*').remove();
      try { d3.select('#hexbin-compare').selectAll('svg').remove(); } catch { }
      this.svg_compare = null as any;
      this.g_compare = null as any;

      // Reset compare legend / scales so reopening starts fresh
      try {
        this.currentLegendDomainCompare = [];
        this.currentCompareLegendType = 'categorical';
        // reset compare ordinal domain
        if (this.colorScaleCompare && typeof this.colorScaleCompare.domain === 'function') {
          this.colorScaleCompare.domain([] as any);
        }
        // reset compare continuous domain
        if (this.continuousColorScaleCompare && typeof this.continuousColorScaleCompare.domain === 'function') {
          this.continuousColorScaleCompare.domain([0, 1]);
        }
      } catch (e) {
        // ignore
      }
    }
  }

  private zoomed(event: d3.D3ZoomEvent<SVGSVGElement, unknown>): void {
    this.currentTransform = event.transform;
    this.g.attr('transform', event.transform.toString());
  }


  private loadAndRenderData(dataPath: string): void {
    // Load GeoJSON data - prepend base URL if it starts with /api/
    const fullUrl = dataPath.startsWith('/api/')
      ? `${this.sessionService.apiUrl}${dataPath}`
      : dataPath;

    // Use fetch instead of d3.json to handle credentials properly
    fetch(fullUrl, { credentials: 'include' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`HTTP Error: ${response.status} ${response.statusText}`);
        }
        return response.json();
      })
      .then((data: GeoJsonData) => {
        console.log('Data loaded from:', fullUrl, data);

        if (!data) {
          throw new Error('Failed to load GeoJSON data');
        }

        this.geoDataService.setData(data);

        // Infer properties for coloring dropdown
        this.features = data.features;

        if (data.meta) {
          // save all features for potential downsampling in Xenium datasets
          this.fullFeatures = data.features;

          this.isXenium =
            (data.meta && data.meta["data_type"] === "xenium") ||
            data.features.length > 50000;
          if (this.svg) {
            if (this.isXenium) {
              (this.svg as any).on('dblclick.zoom', null);
            }
          }
          // set base
          if (this.isXenium) {
            const target = 10000;
            const step = Math.ceil(this.fullFeatures.length / target);
            this.features = this.fullFeatures.filter((_, i) => i % step === 0);

            // window only for xenium
            if (this.isXenium) {

              this.initDetailWindow();
              this.hideDetailWindow();
              this.bindDetailWindowInteractions();
            }

                    } else {
                        this.features = this.fullFeatures;
                    }
                    this.meta = data.meta;
                    this.dgeaReady = !!this.meta?.['dgea']?.['cell_type'];
                    if (this.dgeaReady) {
                      this.initDgeaSelection();
                    }
                    const leidenClusterAnnotations = this.meta?.['leiden_cluster_annotations'];
                    if (leidenClusterAnnotations && typeof leidenClusterAnnotations === 'object') {
                      this.clusterCount = Object.keys(leidenClusterAnnotations).length;
                    }

                    const interval = this.meta?.['interval'];
                    if (Array.isArray(interval) && interval.length > 0) {
                      this.maxInterval = interval.length - 1;
                    }
                    this.selectedRegulatoryScore =
                        this.meta['grn_score_names']?.[0] || null;
                    this.geneSetsGenie3 = this.meta['genie_genesets'] || {};
                    this.geneSetsSponge = this.meta['sponge_genesets'] || {};
                    this.selectedGeneSetGenie3 =
                        Object.keys(this.meta['genie_genesets'] || {})[0] || null;
                    this.selectedGeneSetSponge =
                        Object.keys(this.meta['sponge_genesets'] || {})[0] || null;
                    this.previousGeneSetGenie3 =
                        Object.keys(this.meta['genie_genesets'] || {})[0] || null;
                    this.previousGeneSetSponge =
                        Object.keys(this.meta['sponge_genesets'] || {})[0] || null;
                }

        let firstProps = this.features[0]?.properties || {};

        // start with whatever is available on the first feature, but also make
        // sure that the set contains the "known" columns we care about so that
        // they are present (albeit disabled) until the corresponding data
        // arrives.  This mirrors the behaviour you wanted for regulatory scores
        // without hardcoding a single property.
        const scoreKeys = ['leiden', 'regulatory_scores', 'gene_expression'];
        const lianaKeys = [
          'ligand_receptor_relationships',
          'cell_comp_tf_activity_similarity',
          'tf_activity',
          'pathway_activity',
        ];

        const allKeys = new Set<string>(
          Object.keys(firstProps).filter((k) => {
            const val = firstProps[k];
            return typeof val === 'string' || typeof val === 'number';
          })
        );
        scoreKeys.forEach((k) => allKeys.add(k));
        lianaKeys.forEach((k) => allKeys.add(k));

        this.leidenCentralityProps.forEach((k) => allKeys.add(k));
        this.colorableProperties = Array.from(allKeys);
        console.log('Initial colourable properties (including potential ones):',
          this.colorableProperties);

        // Alphabetical order
        this.colorableProperties.sort((a, b) => a.localeCompare(b));

        // Group similar properties together
        const chromvarKeys = ['chromvar_total_sum'];
        this.groupedProperties = [
          { key: 'Scores', value: this.colorableProperties.filter((p) => scoreKeys.includes(p)) },
          { key: 'LIANA+', value: this.colorableProperties.filter((p) => lianaKeys.includes(p)) },
          { key: 'ChromVAR', value: this.colorableProperties.filter((p) => chromvarKeys.includes(p)) },
          {
            key: 'Other', value: this.colorableProperties.filter(
              (p) => !scoreKeys.includes(p) && !lianaKeys.includes(p) && !chromvarKeys.includes(p)
            )
          },
        ];

        // compute availability for everything we've decided to show
        const candidates = new Set(this.colorableProperties);
        // ensure some non-colorable but still important keys are included
        candidates.add('regulatory_scores');
        this.propertyAvailability = {};
        candidates.forEach((prop) => {
          if (prop === 'regulatory_scores') {
            this.propertyAvailability[prop] = this.hasRegulatoryScoresData();
          }else if (this.leidenCentralityProps.includes(prop)) {
            this.propertyAvailability[prop] = this.features.some((f) => {
              const val = this.getLeidenClusterAnnotation(f.properties.leiden)?.centrality?.[prop];
              return val !== undefined && val !== null && val !== '';
            });
          } else {
            this.propertyAvailability[prop] = this.features.some((f) => {
              const val = f.properties ? f.properties[prop] : undefined;
              return val !== undefined && val !== null && val !== '';
            });
          }
        });
        console.log('Property availability map:', this.propertyAvailability);

        // remove unavailable properties from the dropdown entirely so the user
        // cannot select them and they don't appear as the default
        this.colorableProperties = this.colorableProperties.filter((p) =>
          this.propertyAvailability[p]
        );
        console.log('Filtered colourableProperties:', this.colorableProperties);

        // rebuild grouped properties after filtering to exclude unavailable items
        this.groupedProperties = [
          { key: 'Scores', value: this.colorableProperties.filter((p) => scoreKeys.includes(p)) },
          { key: 'LIANA+', value: this.colorableProperties.filter((p) => lianaKeys.includes(p)) },
          {
            key: 'Other', value: this.colorableProperties.filter(
              (p) => !scoreKeys.includes(p) && !lianaKeys.includes(p)
            ),
          },
        ];

        if (this.propertyAvailable('regulatory_scores')) {
          this.colorByProperty = 'regulatory_scores';
        } else if (this.propertyAvailable('cell_type')) {
          this.colorByProperty = 'cell_type';
        } else {
          // just pick the first remaining option (already filtered to available
          // ones).  any value here is acceptable as 'random' default.
          this.colorByProperty = this.colorableProperties[0] || '';
        }

        this.currentLegendType = this.isContinuousScale() ? 'continuous' : 'categorical';

        const width = 1200;
        const height = 1000;

        const projection = d3.geoIdentity().fitSize([width, height], {
          type: 'FeatureCollection',
          features: (this.isXenium ? this.fullFeatures : this.features),
        });

        // Create a geoPath generator with the projection
        this.currentPathGenerator = d3.geoPath<CellFeature>().projection(projection);
        const pathGenerator = this.currentPathGenerator;

        // Precompute projected centroids for Xenium detail window (performance)
        if (this.isXenium) {
          this.fullFeatures.forEach((f) => {
            if (!(f.properties as any).__centroidProjected) {
              const c = this.currentPathGenerator!.centroid(f as any);
              (f.properties as any).__centroidProjected = c;
            }
          });
        }



        // Ziel-Layer auswählen (Xenium = baseLayer, Visium = g)
        if (this.isXenium) {
          this.baseLayer
              .style('cursor', 'default')
              .style('pointer-events', 'none')
              .selectAll<SVGPathElement, CellFeature>('path')
              .data(this.features, (d: any) => d.properties.barcode)
              .join('path')
              .attr('d', (d: CellFeature) => pathGenerator(d) || '')
              .attr('fill', (d: CellFeature) => {
                  const value = this.leidenCentralityProps.includes(this.colorByProperty)
                    ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[this.colorByProperty]
                    : d.properties?.[this.colorByProperty];
                  if (this.currentLegendType === 'categorical') {
                      return this.colorScale(String(value));
                  } else {
                      const num = this.toNumber(value);
                      return Number.isFinite(num)
                          ? this.continuousColorScale(num)
                          : '#ccc';
                  }
              })
              .style('opacity', 0.8);
        } else {
          this.g
              .style('cursor', 'pointer')
              .style('pointer-events', null)
              .selectAll<SVGPathElement, CellFeature>('path')
              .data(this.features, (d: any) => d.properties.barcode)
              .join('path')
              .attr('d', (d: CellFeature) => pathGenerator(d) || '')
              .attr('fill', (d: CellFeature) => {
                  const value = this.leidenCentralityProps.includes(this.colorByProperty)
                    ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[this.colorByProperty]
                    : d.properties?.[this.colorByProperty];
                  if (this.currentLegendType === 'categorical') {
                      return this.colorScale(String(value));
                  } else {
                      const num = this.toNumber(value);
                      return Number.isFinite(num)
                          ? this.continuousColorScale(num)
                          : '#ccc';
                  }
              })
              .style('opacity', 0.8)
              .on('mouseover', (event, d) => this.mouseOver(event, d))
              .on('mouseleave', (event, d) => this.mouseLeave(event, d))
              .on('click', (event, d) => this.openSidenav(event, d));
        }

        this.onColorbyPropertyChange();
        setTimeout(() => {
          this.isLoadingHexagons = false;
        }, 0);
      })
      .catch((error) => {
        console.error('Error loading or rendering data:', error);
        this.isLoadingHexagons = false;
      });
  }


  private toNumber(v: unknown): number {
    if (typeof v === 'number') return v;
    if (typeof v === 'string') {
      const n = parseFloat(v);
      return Number.isFinite(n) ? n : NaN;
    }
    return NaN;
  }

  public onColorbyPropertyChange(): void {
    console.log('[onColorbyPropertyChange] colorByProperty changed to:', this.colorByProperty);

    if (this.colorByProperty === 'regulatory_scores') {
      if (
        this.selectedRegulatoryScore?.endsWith('genie3') &&
        this.selectedGeneSetGenie3
      ) {
        console.log(this.selectedRegulatoryScore, this.selectedGeneSetGenie3);
        this.fetchAndUpdate(
          this.selectedRegulatoryScore,
          this.selectedGeneSetGenie3,
        );
        this.updateSubgraphGenie3();
        this.updateSubgraphSponge();

      } else if (
        this.selectedRegulatoryScore?.endsWith('sponge') &&
        this.selectedGeneSetSponge
      ) {
        this.fetchAndUpdate(
          this.selectedRegulatoryScore,
          this.selectedGeneSetSponge,
        );
        this.updateSubgraphGenie3();
        this.updateSubgraphSponge();
      }
    }
    if (this.compareMode) {
      console.log('[onColorbyPropertyChange] Updating compare hexagons');
      this.updateHexColors('#hexbin-compare');
    }
    console.log('[onColorbyPropertyChange] Updating main hexagons');
    this.updateHexColors();
  }

  /**
   * Determine whether the given view/property should be treated as a continuous numeric scale.
   * If `features` is provided, that dataset is used for the test (used for compare view). Otherwise
   * the main `this.features` is used.
   */
  isContinuousScale(view?: string, features?: CellFeature[]) {
    if (!view) {
      view = this.colorByProperty;
    }
    const sourceFeatures = Array.isArray(features) ? features : this.features || [];

        const valuesRaw = sourceFeatures.map((f) => {
            if (this.leidenCentralityProps.includes(view)) {
                const clusterAnnotation = this.getLeidenClusterAnnotation(f.properties.leiden);
                return clusterAnnotation?.centrality?.[view];
            }
            return f.properties[view];;
        });

    const numericValues = valuesRaw.map((v) => this.toNumber(v));
    const allNumbers = numericValues.every((n) => Number.isFinite(n));

    // Check if all values are integers (for categorical treatment)
    const allIntegers = allNumbers && numericValues.every((n) => Number.isInteger(n));

    // Check if we have a reasonable number of unique integer values for categorical treatment (here 20)
    const uniqueIntegerCount = allIntegers ? new Set(numericValues).size : 0;
    const shouldTreatAsCategorical = allIntegers && uniqueIntegerCount <= 20;

    return allNumbers && !shouldTreatAsCategorical && numericValues.length > 0;
  }

  private getViewVariablesToUpdate(containerName: string) {
    const isMainView = containerName === '#hexbin';
    const viewToUpdate = isMainView ? this.colorByProperty : this.selectedCompareView;
    const featuresToUpdate = isMainView ? this.features : (this.dataCompare?.features || []);
    const gToUpdate = isMainView ? this.g : this.g_compare;
    const gPathsToUpdate = isMainView ? this.g_paths : this.g_paths_compare;
    const ordinalScaleToUpdate = isMainView ? this.colorScale : this.colorScaleCompare;
    const continuousScaleToUpdate = isMainView ? this.continuousColorScale : this.continuousColorScaleCompare;
    const isContinuous = this.isContinuousScale(viewToUpdate, featuresToUpdate);
    const legendContainerName = isMainView ? 'svg-legend' : 'svg-legend-compare';
    const legendGradientName = isMainView ? 'svg-legend-gradient' : 'svg-legend-gradient-compare';
    return {
      isMainView,
      view: viewToUpdate,
      svg: isMainView ? this.svg : this.svg_compare,
      features: featuresToUpdate,
      g: gToUpdate,
      g_paths: gPathsToUpdate,
      ordinal: ordinalScaleToUpdate,
      continuous: continuousScaleToUpdate,
      isContinuous,
      legendContainerName,
      legendGradientName,
      getLegendDomain: () => isMainView ? this.currentLegendDomain : this.currentLegendDomainCompare,
      setLegendDomain: (v: any[]) => { if (isMainView) this.currentLegendDomain = v; else this.currentLegendDomainCompare = v; },
      getLegendType: () => isMainView ? this.currentLegendType : this.currentCompareLegendType,
      setLegendType: (t: 'continuous' | 'categorical') => { if (isMainView) this.currentLegendType = t; else this.currentCompareLegendType = t; }
    } as const;

  }

  public updateHexColors(containerName?: string): void {

    if (!containerName) {
      containerName = '#hexbin';
    }

      const viewVariablesToUpdate = this.getViewVariablesToUpdate(containerName);
      let valuesRaw;

      console.log('[updateHexColors] container=', containerName, 'view=', viewVariablesToUpdate.view, 'featuresForTest=', viewVariablesToUpdate.features.length || 0, 'isContinuous=', viewVariablesToUpdate.isContinuous);

      // Debug: check if the view property exists in features
      if (viewVariablesToUpdate.features.length > 0) {
        const firstFeature = viewVariablesToUpdate.features[0];
        const hasProperty = this.leidenCentralityProps.includes(viewVariablesToUpdate.view)
          ? true
          : viewVariablesToUpdate.view in (firstFeature.properties || {});
        console.log('[updateHexColors] First feature properties keys:', Object.keys(firstFeature.properties || {}));
        console.log('[updateHexColors] Looking for property:', viewVariablesToUpdate.view, '- Exists:', hasProperty);
        if (!hasProperty) {
          console.warn('[updateHexColors] Property not found in features! Available:', Object.keys(firstFeature.properties || {}));
        }
      }

    this.resetClusterExtension();

    if (this.selectedCell && this.selectedCluster) {
      this.selectedCluster = null;
      this.clusterCells = [];
      this.clusterCellTypes = [];
      this.clusterCentralityAvg = {
        degree_centrality: 0,
        average_clustering: 0,
        closeness_centrality: 0,
      };
    }
    if (this.selectedCell) this.selectedCell = null;


    const layerToColor = (
      (containerName === '#hexbin' && this.isXenium)
        ? this.baseLayer
        : viewVariablesToUpdate.g
    ) as unknown as d3.Selection<SVGGElement, any, any, any>;

    const sel = layerToColor
      .selectAll('path')
      .data(viewVariablesToUpdate.features);

        valuesRaw = viewVariablesToUpdate.features.map((f) => {
            if (this.leidenCentralityProps.includes(viewVariablesToUpdate.view)) {
                return this.getLeidenClusterAnnotation(f.properties.leiden)?.centrality?.[viewVariablesToUpdate.view];
            }
            return f.properties[viewVariablesToUpdate.view];
        });

    const numericValues = valuesRaw.map((v) => this.toNumber(v));

    if (viewVariablesToUpdate.isContinuous) {
      // continuous scale - only if not integers or too many unique integers
      let min = Math.min(...numericValues);
      let max = Math.max(...numericValues);
      if (min === max) {
        const eps = min === 0 ? 1 : Math.abs(min) * 0.01;
        min -= eps;
        max += eps;
      }


      viewVariablesToUpdate.continuous.domain([min, max]);

      viewVariablesToUpdate.setLegendDomain([min, max]);
      viewVariablesToUpdate.setLegendType('continuous');

            sel
                .transition()
                .duration(300)
                .attr('stroke-width', 1)
                .attr('stroke', 'transparent')
                .attr('fill', (d) => {
                    const raw = this.leidenCentralityProps.includes(viewVariablesToUpdate.view)
                        ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[viewVariablesToUpdate.view]
                        : d.properties[viewVariablesToUpdate.view];
                    const n = this.toNumber(raw);
                    return Number.isFinite(n)
                        ? viewVariablesToUpdate.continuous(n)
                        : '#ccc';

        });
    } else {
      // categorical scale - for non-numeric, integers with few unique values, or mixed data
      const domain = [...new Set(valuesRaw.map((v: any) => String(v)))];

      viewVariablesToUpdate.ordinal.domain(domain);
      viewVariablesToUpdate.setLegendDomain(domain);
      viewVariablesToUpdate.setLegendType('categorical');

            sel
                .transition()
                .duration(300)
                .attr('stroke-width', 1)
                .attr('stroke', 'transparent')
                .attr('fill', (d) => {
                    const raw = this.leidenCentralityProps.includes(viewVariablesToUpdate.view)
                        ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[viewVariablesToUpdate.view]
                        : d.properties[viewVariablesToUpdate.view];
                    return viewVariablesToUpdate.ordinal(String(raw));
                });
        }
        if (
            containerName === '#hexbin' &&
            this.isXenium &&
            this.detailVisible &&
            this.detailScreenPos
        ) {
            this.showDetailWindowAt(this.detailScreenPos.x, this.detailScreenPos.y);
            this.updateDetailAtScreenPos(this.detailScreenPos.x, this.detailScreenPos.y);
        }
        this.renderLegend(containerName);
    }

  public updateSubgraphGenie3(): void {
    const token = this.nextRequestToken('genie3');
    console.log('Updating AUCELL graph for Genie3...');
    this.isLoadingGenie3 = true;
    d3.select('#aucell_graph_genie3').selectAll('*').remove();

    if (!this.selectedGeneSetGenie3 || !this.genie3Network) {
      this.isLoadingGenie3 = false;
      return;
    }

    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/geneset_connections_genie?gene_set_name=${encodeURIComponent(this.selectedGeneSetGenie3 ? this.selectedGeneSetGenie3 : '')}`,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          if (token !== this.requestTokens['genie3']) return;
          {
            const payload: any = res;
            const data = payload['connections'] as { regulatoryGene: string; targetGene: string; weight: number }[];
            this.genie3Network = data.map((d) => ({
              source: d.regulatoryGene,
              target: d.targetGene,
              weight: d.weight,
            }));

            // slider_data may be returned as an object or an array; handle both safely
            const sliderData: any = payload['slider_data'];

            if (sliderData && typeof sliderData === 'object') {
              this.genie3SliderData = {
                step: sliderData.step || 1,
                min_border: sliderData.min_border || 0,
                max_border: sliderData.max_border || 100,
                default_value: sliderData.default_value || 50,
              };
              this.genie3WeightCutoff = this.genie3SliderData.default_value;
            }

            this.visualizeGenie3Subgraph();

          }
        },
        error: (err) => {
          if (this.requestTokens['sponge'] !== this.requestTokens['sponge']) {
            console.error(
              `[Backend] Failed to load Sponge Connections for["${this.selectedGeneSetSponge}]`,
              err,
            );
          }
        },
      });
  }

  public visualizeGenie3Subgraph(): void {

    this.isLoadingGenie3 = true;
    d3.select('#aucell_graph_genie3').html('');

    if (!this.selectedGeneSetGenie3 || !this.genie3Network) {
      this.isLoadingGenie3 = false;
      return;
    }

    let regulator = this.selectedGeneSetGenie3;
    let targets = this.geneSetsGenie3[regulator] || [];

    let nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    let edges: { source: string; target: string; weight: number }[] = [];

    let candidateEdges: { source: string; target: string; weight: number }[] = [];
    let slicedEdges: { source: string; target: string; weight: number }[] = [];


    candidateEdges = this.genie3Network.filter((edge) => edge.weight > this.genie3WeightCutoff).map((e) => ({
      source: String(e.source),
      target: String(e.target),
      weight: e.weight,
    }));

    candidateEdges.sort((a, b) => b.weight - a.weight);
    slicedEdges = candidateEdges.slice(0, this.genie3MinEdges);

    // Infer nodes from edges
    const nodeSet = new Set<string>();

    slicedEdges.forEach((edge) => {
      if (!(nodeSet.has(edge.source) && nodeSet.has(edge.target))) {
        console.log("Adding edge nodes:", edge.source, edge.target);
        nodeSet.add(edge.source);
        nodeSet.add(edge.target);
      }
    });

    // Add regulator to node set
    nodeSet.add(regulator);

    // Readd edges that passed the cutoff
    candidateEdges.forEach((edge) => {
      if (nodeSet.has(edge.source) && nodeSet.has(edge.target)) {
        edges.push(edge);
      }
    });

    // Add nodes with groups
    nodeSet.forEach((nodeId) => {
      if (nodeId === regulator) {
        nodes.push({ id: nodeId, group: 0 });
      } else if (targets.includes(nodeId)) {
        nodes.push({ id: nodeId, group: 1 });
      } else {
        nodes.push({ id: nodeId, group: 2 });
      }
    });

    console.log(nodes);

    edges.push(...candidateEdges);

    // Create the graph
    let graph = {
      nodes: nodes.filter((node) => node.id && node.id.length > 0),
      edges: edges.filter(
        (edge) =>
          nodes.some((node) => node.id === edge.source) &&
          nodes.some((node) => node.id === edge.target),
      ),
    };

    console.log(graph)

    // Create the graph visualization
    const width = this.genie3Width;
    const height = 300;

    const svg = d3
      .select('#aucell_graph_genie3')
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .style('background-color', '#f8f9fa');

    // Draw links (edges)
    const link = svg
      .append('g')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.6)
      .selectAll('line')
      .data(graph.edges)
      .enter()
      .append('line')
      .attr('stroke-width', (d: any) =>
        Math.max(1, Math.sqrt(d.weight) * 10),
      )
      .attr('stroke', (d: any) => {
        // Color edges based on weight
        const intensity = Math.min(d.weight * 10, 1);
        return d3.interpolateReds(0.3 + intensity * 0.7);
      });

    // Draw nodes
    const node = svg
      .append('g')
      .attr('stroke', '#fff')
      .attr('stroke-width', 1.5)
      .selectAll('circle')
      .data(graph.nodes)
      .enter()
      .append('circle')
      .attr('r', (d: any) => {
        switch (d.group) {
          case 0:
            return 15; // regulator
          case 1:
            return 12; // targets
          case 2:
            return 8; // neighbors
          case 3:
            return 6;
          default:
            return 10;
        }
      })
      .attr('fill', (d: any) => {
        switch (d.group) {
          case 0:
            return '#e41a1c'; // regulator - red
          case 1:
            return '#377eb8'; // targets - blue
          case 2:
            return '#4daf4a'; // neighbors - green
          case 3:
            return '#ff7f00'; // high-weight - orange
          default:
            return '#999';
        }
      });

    // Add labels
    const labels = svg
      .append('g')
      .selectAll('text')
      .data(graph.nodes)
      .enter()
      .append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '.35em')
      .style('font-size', (d: any) => (d.group === 0 ? '12px' : '10px'))
      .style('font-weight', (d: any) =>
        d.group === 0 ? 'bold' : 'normal',
      )
      .style('fill', '#333')
      .text((d: any) =>
        d.id.length > 8 ? d.id.substring(0, 8) + '...' : d.id,
      );

    // Initialize simulation with stronger forces
    const simulation = d3
      .forceSimulation(graph.nodes)
      .force(
        'link',
        d3
          .forceLink(graph.edges)
          .id((d: any) => d.id)
          .distance(30)
          .strength(0.5),
      )
      .force('charge', d3.forceManyBody().strength(-500))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(30))
      .force('boundary', () => {
        graph.nodes.forEach((node: any) => {
          const radius =
            node.group === 0
              ? 15
              : node.group === 1
                ? 12
                : node.group === 2
                  ? 8
                  : 6;

          node.x = Math.max(radius, Math.min(width - radius, node.x));

          node.y = Math.max(radius, Math.min(height - radius, node.y));
        });
      });

    simulation.on('tick', () => {
      link
        .attr('x1', (d: any) => d.source.x)
        .attr('y1', (d: any) => d.source.y)
        .attr('x2', (d: any) => d.target.x)
        .attr('y2', (d: any) => d.target.y);

      node.attr('cx', (d: any) => d.x).attr('cy', (d: any) => d.y);

      labels.attr('x', (d: any) => d.x).attr('y', (d: any) => d.y);
    });

    this.isLoadingGenie3 = false;


  }


  public updateSubgraphSponge(): void {
    const token = this.nextRequestToken('sponge');
    console.log('Updating AUCELL graph for Sponge...');
    d3.select('#aucell_graph_sponge').selectAll('*').remove();

    if (!this.selectedGeneSetSponge || !this.spongeNetwork) {
      return;
    }

    this.isLoadingSponge = true;
    const regulator = this.selectedGeneSetSponge;

    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/geneset_connections_sponge?gene_set_name=${encodeURIComponent(this.selectedGeneSetSponge ? this.selectedGeneSetSponge : '')}`,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          if (token !== this.requestTokens['sponge']) return;
          {
            const payload: any = res;

            const data = payload['connections'] as {
              geneA: string;
              geneB: string;
              'p.adj': number;
              mscor: number;
            }[];
            const sliderData: any = payload['slider_data'];

            if (sliderData && typeof sliderData === 'object') {
              this.spongeSliderData = {
                step: sliderData.step || 0.01,
                min_border: sliderData.min_border || 0,
                max_border: sliderData.max_border || 100,
                default_value: sliderData.default_value || 50,
              };
            }

            this.isLoadingSponge = true;

            console.log('Sponge Network:', data);

            this.spongeNetwork = data.map((d) => ({
              source: d.geneA,
              target: d.geneB,
              p_adjusted: d['p.adj'],
              mscore: d['mscor'],
            }));

            this.visualizeSpongeSubgraph();
          }
        },
        error: (err) => {
          if (this.requestTokens['sponge'] !== this.requestTokens['sponge']) {
            console.error(
              `[Backend] Failed to load Sponge Connections for["${this.selectedGeneSetSponge}]`,
              err,
            );
          }
        },
      });
  }

  public visualizeSpongeSubgraph(): void {

    this.isLoadingSponge = true;
    d3.select('#aucell_graph_sponge').html('');

    if (!this.selectedGeneSetSponge || !this.spongeNetwork) {
      this.isLoadingSponge = false;
      return;
    }

    let regulator = this.selectedGeneSetSponge;
    let targets = this.geneSetsSponge[regulator] || [];

    let nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    let edges: { source: string; target: string; p_adjusted: number }[] = [];

    let candidateEdges: { source: string; target: string; p_adjusted: number }[] = [];
    let slicedEdges: { source: string; target: string; p_adjusted: number }[] = [];


    candidateEdges = this.spongeNetwork.filter((edge) => edge.p_adjusted < this.spongePValueCutoff).map((e) => ({
      source: String(e.source),
      target: String(e.target),
      p_adjusted: e.p_adjusted,
    }));

    candidateEdges.sort((a, b) => a.p_adjusted - b.p_adjusted);
    slicedEdges = candidateEdges.slice(0, this.spongeMinEdges);

    // Infer nodes from edges
    const nodeSet = new Set<string>();

    slicedEdges.forEach((edge) => {
      if (!(nodeSet.has(edge.source) && nodeSet.has(edge.target))) {
        console.log("Adding edge nodes:", edge.source, edge.target);
        nodeSet.add(edge.source);
        nodeSet.add(edge.target);
      }
    });

    // Add regulator to node set
    nodeSet.add(regulator);

    // Add nodes with groups
    nodeSet.forEach((nodeId) => {
      if (nodeId === regulator) {
        nodes.push({ id: nodeId, group: 0 });
      } else if (targets.includes(nodeId)) {
        nodes.push({ id: nodeId, group: 1 });
      } else {
        nodes.push({ id: nodeId, group: 2 });
      }
    });

    // Readd edges of the nodes that passed the cutoff
    candidateEdges.forEach((edge) => {
      if (nodeSet.has(edge.source) && nodeSet.has(edge.target)) {
        edges.push(edge);
      }
    });

    edges.push(...candidateEdges);

    // Create the graph
    let graph = {
      nodes: nodes.filter((node) => node.id && node.id.length > 0),
      edges: edges.filter(
        (edge) =>
          nodes.some((node) => node.id === edge.source) &&
          nodes.some((node) => node.id === edge.target),
      ),
    };

    // Create the graph visualization
    const width = this.spongeWidth;
    const height = 300;

    const svg = d3
      .select('#aucell_graph_sponge')
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .style('background-color', '#f8f9fa');

    // Draw links (edges)
    const link = svg
      .append('g')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.6)
      .selectAll('line')
      .data(graph.edges)
      .enter()
      .append('line')
      .attr('stroke-width', (d: any) =>
        Math.max(1, Math.sqrt(d.p_adjusted) * 10),
      )
      .attr('stroke', (d: any) => {
        // Color edges based on p_adjusted
        const intensity = Math.min(d.p_adjusted * 10, 1);
        return d3.interpolateReds(0.3 + intensity * 0.7);
      });

    // Draw nodes
    const node = svg
      .append('g')
      .attr('stroke', '#fff')
      .attr('stroke-width', 1.5)
      .selectAll('circle')
      .data(graph.nodes)
      .enter()
      .append('circle')
      .attr('r', (d: any) => {
        switch (d.group) {
          case 0:
            return 15; // regulator
          case 1:
            return 12; // targets
          case 2:
            return 8; // neighbors
          case 3:
            return 6;
          default:
            return 10;
        }
      })
      .attr('fill', (d: any) => {
        switch (d.group) {
          case 0:
            return '#e41a1c'; // regulator - red
          case 1:
            return '#377eb8'; // targets - blue
          case 2:
            return '#4daf4a'; // neighbors - green
          case 3:
            return '#ff7f00'; // high-weight - orange
          default:
            return '#999';
        }
      });

    // Add labels
    const labels = svg
      .append('g')
      .selectAll('text')
      .data(graph.nodes)
      .enter()
      .append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '.35em')
      .style('font-size', (d: any) => (d.group === 0 ? '12px' : '10px'))
      .style('font-weight', (d: any) =>
        d.group === 0 ? 'bold' : 'normal',
      )
      .style('fill', '#333')
      .text((d: any) =>
        d.id.length > 8 ? d.id.substring(0, 8) + '...' : d.id,
      );

    // Initialize simulation with stronger forces
    const simulation = d3
      .forceSimulation(graph.nodes)
      .force(
        'link',
        d3
          .forceLink(graph.edges)
          .id((d: any) => d.id)
          .distance(20)
          .strength(0.5),
      )
      .force('charge', d3.forceManyBody().strength(-500))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(30))
      .force('boundary', () => {
        graph.nodes.forEach((node: any) => {
          const radius =
            node.group === 0
              ? 15
              : node.group === 1
                ? 12
                : node.group === 2
                  ? 8
                  : 6;

          node.x = Math.max(radius, Math.min(width - radius, node.x));

          node.y = Math.max(radius, Math.min(height - radius, node.y));
        });
      });

    simulation.on('tick', () => {
      link
        .attr('x1', (d: any) => d.source.x)
        .attr('y1', (d: any) => d.source.y)
        .attr('x2', (d: any) => d.target.x)
        .attr('y2', (d: any) => d.target.y);

      node.attr('cx', (d: any) => d.x).attr('cy', (d: any) => d.y);

      labels.attr('x', (d: any) => d.x).attr('y', (d: any) => d.y);
    });

    this.isLoadingSponge = false;

    console.log('Network visualization complete');
  }


  public onSliderChangeSponge(): void {
    if (this.selectedGeneSetSponge) {
      this.isLoadingSponge = true;
      setTimeout(() => this.visualizeSpongeSubgraph(), 50);
    }
  }


  public onSliderChangeGenie3(): void {
    if (this.selectedGeneSetGenie3) {
      this.isLoadingGenie3 = true;
      setTimeout(() => this.visualizeGenie3Subgraph(), 50);
    }
  }

  public analyzeGeneSetInGProfiler(): void {
    if (!this.selectedGeneSetGenie3 || !this.geneSetsGenie3) {
      console.warn('No Genie3 gene set selected');
      return;
    }

    const regulator = this.selectedGeneSetGenie3;
    const targets = this.geneSetsGenie3[regulator] || [];

    const allGenes = [regulator, ...targets];

    console.log('Analyzing all Genie3 genes in gProfiler:', allGenes);

    const gProfilerUrl = this.generateGProfilerUrl(allGenes);

    if (gProfilerUrl) {
      // Open in new tab/window
      window.open(gProfilerUrl, '_blank');
    } else {
      console.warn(
        'Could not generate gProfiler URL for Genie3 gene set:',
        this.selectedGeneSetGenie3,
      );
    }
  }

  // Make sure you also have the generateGProfilerUrl method
  private generateGProfilerUrl(geneIds: string[]): string | null {
    if (!geneIds || geneIds.length === 0) {
      return null;
    }

    // Join gene IDs with newlines (gProfiler expects one gene per line)
    const geneList = geneIds.join('\n');

    // Base gProfiler URL for functional enrichment analysis
    const baseUrl = 'https://biit.cs.ut.ee/gprofiler/gost';

    // URL encode the gene list
    const encodedGenes = encodeURIComponent(geneList);

    // Construct the full URL with parameters including auto-run
    const gProfilerUrl = `${baseUrl}?organism=hsapiens&query=${encodedGenes}&sources=GO:MF,GO:BP,GO:CC,KEGG,REAC&user_threshold=0.05&significance_threshold_method=fdr&ordered=false&exclude_iea=false&measure_underrepresentation=false&evcodes=false&domain_scope=annotated&numeric_ns=ENTREZGENE_ACC&background=&run_query=1`;

    return gProfilerUrl;
  }

  private mouseOver(event: MouseEvent, d: CellFeature): void {
    d3.selectAll('.Country')
      .transition()
      .duration(200)
      .style('opacity', 0.5)
      .attr('stroke', 'transparent');

    d3.select(event.target as SVGElement)
      .transition()
      .duration(200)
      .style('opacity', 0.8)
      .attr('stroke', 'black');
  }

  private mouseLeave(event: MouseEvent, d: CellFeature): void {
    if (
      this.selectedCell &&
      (d.properties.barcode === this.selectedCell.properties.barcode ||
        (this.colorByProperty === 'leiden' &&
          d.properties.leiden === this.selectedCell.properties.leiden))
    )
      return;
    if (
      this.selectedCellCompare &&
      (d.properties.barcode === this.selectedCellCompare.properties.barcode ||
        (this.selectedCompareView === 'leiden' &&
          d.properties.leiden === this.selectedCellCompare.properties.leiden))
    )
      return;
    d3.selectAll('.Country')
      .transition()
      .duration(200)
      .style('opacity', 0.8)
      .attr('stroke', 'transparent');

    d3.select(event.target as SVGElement)
      .transition()
      .duration(200)
      .attr('stroke', 'transparent');
  }

  public openSidenav(event: MouseEvent, cell: CellFeature): void {
    this.resetClusterExtension();
    this.selectedCell = cell;
    if (this.colorByProperty === 'regulatory_scores') {
      this.getRegulatoryScoresforSpots(cell.properties.barcode)
    }
    if (this.colorByProperty === 'leiden') {
      this.openClusterSidenav(cell.properties.leiden);
      this.extendCluster(cell.properties.leiden);
    } else {
      d3.select(event.target as SVGElement)
        .transition()
        .attr('stroke', 'black');
    }

    // Neighborhood enrichment will render when the tab is viewed
    // via onTabChange handler
    setTimeout(() => this.updateSubgraphGenie3(), 0);
  }

  public openSidenavCompare(event: MouseEvent, cell: CellFeature): void {
    this.selectedCellCompare = cell;
    if (this.selectedCompareView === 'regulatory_scores') {
      this.getRegulatoryScoresforSpots(cell.properties.barcode)
    }
    d3.select(event.target as SVGElement)
      .transition()
      .attr('stroke', 'black');
  }

  public openClusterSidenav(clusterId: number): void {
    this.selectedCluster = clusterId;
    this.clusterCells = this.features.filter(
      (cell) => cell.properties.leiden === clusterId,
    );
    this.calculateClusterStats();

    // Initialize co-occurrence table for this cluster
    this.updateCoOccurrenceTable();

    if (this.clusterCells.length > 0) {
      this.selectedCell = this.clusterCells[0];
      // Neighborhood enrichment will render when the Cluster Information tab is viewed
      // via onTabChange handler
      setTimeout(() => this.updateSubgraphGenie3(), 100);
      setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
      // setTimeout(() => this.renderFootprintPlots(), 0);
    }
  }

  public selectCluster(clusterId: number): void {
    this.selectedCluster = clusterId;
    this.clusterCells = this.features.filter(
      (cell) => cell.properties.leiden === clusterId,
    );
    this.calculateClusterStats();
    this.updateCoOccurrenceTable();

    if (this.clusterCells.length > 0) {
      this.selectedCell = this.clusterCells[0];
      // Neighborhood enrichment will render when the Cluster Information tab is viewed
      // via onTabChange handler
      setTimeout(() => this.updateSubgraphGenie3(), 100);
      setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
      // setTimeout(() => this.renderFootprintPlots(), 0);
    }
  }

  public onGeneSetChange(): void {
    if (this.selectedGeneSetGenie3 !== this.previousGeneSetGenie3) {
      this.previousGeneSetGenie3 = this.selectedGeneSetGenie3;

      if (this.selectedGeneSetGenie3) {
        d3.select('#aucell_graph_genie3').selectAll('*').remove();
        this.isLoadingGenie3 = true;
        setTimeout(() => {
          this.updateSubgraphGenie3();
          if (
            this.selectedRegulatoryScore?.endsWith('genie3') &&
            this.selectedGeneSetGenie3
          ) {
            this.fetchAndUpdate(
              this.selectedRegulatoryScore,
              this.selectedGeneSetGenie3,
            );
          }
        }, 100);
      } else {
        // Clear Genie3 graph if no gene set is selected
        d3.select('#aucell_graph_genie3').selectAll('*').remove();
        this.isLoadingGenie3 = false;
      }
    }

    if (this.selectedGeneSetSponge !== this.previousGeneSetSponge) {
      this.previousGeneSetSponge = this.selectedGeneSetSponge;
      if (this.selectedGeneSetSponge) {
        d3.select('#aucell_graph_sponge').selectAll('*').remove();
        console.log('Updating Sponge graph for:', this.selectedGeneSetSponge);
        console.log(
          'Sponge targets available:',
          this.geneSetsSponge[this.selectedGeneSetSponge]?.length || 0,
        );
        this.isLoadingSponge = true;
        setTimeout(() => {
          this.updateSubgraphSponge();
          if (
            this.selectedRegulatoryScore?.endsWith('sponge') &&
            this.selectedGeneSetSponge
          ) {
            this.fetchAndUpdate(
              this.selectedRegulatoryScore,
              this.selectedGeneSetSponge,
            );
          }
        }, 100);
      } else {
        // Clear Sponge graph if no gene set is selected
        d3.select('#aucell_graph_sponge').selectAll('*').remove();
        this.isLoadingSponge = false;
      }
    }
  }

  public selectCellFromCluster(cell: CellFeature): void {
    this.selectedCell = cell;
    // Neighborhood enrichment will render when the Cluster Information tab is viewed
    // via onTabChange handler
    setTimeout(() => this.updateSubgraphGenie3(), 0);
  }

  public closeClusterSidenav(): void {
    this.selectedCluster = null;
    this.clusterCells = [];
    this.clusterCellTypes = [];
    this.coOccurrenceData = []; // Clear co-occurrence data
    this.resetClusterExtension();
    this.updateHexColors();
  }

  private calculateClusterStats(): void {
    if (this.clusterCells.length === 0) return;

    // Calculate cell type distribution using existing cell_type property
    const cellTypeMap = new Map<string, number>();
    this.clusterCells.forEach((cell) => {
      const cellType = cell.properties.cell_type;
      cellTypeMap.set(cellType, (cellTypeMap.get(cellType) || 0) + 1);
    });

    this.clusterCellTypes = Array.from(cellTypeMap.entries())
      .map(([type, count]) => ({
        type,
        count,
        percentage: ((count / this.clusterCells.length) * 100).toFixed(1),
      }))
      .sort((a, b) => b.count - a.count);

        const clusterAnnotation = this.getLeidenClusterAnnotation(this.selectedCluster!);

        if (this.clusterCells.length > 0) {
            this.clusterCentralityAvg = {
                degree_centrality: clusterAnnotation?.centrality?.['degree_centrality'] ?? 0,
                average_clustering: clusterAnnotation?.centrality?.['average_clustering'] ?? 0,
                closeness_centrality: clusterAnnotation?.centrality?.['closeness_centrality'] ?? 0,
            };
        }
    }

  private extendCluster(selectedCluster: number): void {
    this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => {
        // Extending the hexagon size by 1.1 is barely noticeable,
        // plus it's infinitely annoying resetting size when switching
        // away from leiden clustering.
        //if (d.properties.leiden === selectedCluster) {
        //  // Scale the hexagon coordinates outward
        //  return this.getScaledPath(d, 1.1); // 10% larger
        //}
        // Return original path for non-selected hexagons
        const projection = d3.geoIdentity().fitSize([1200, 1000], {
          type: 'FeatureCollection',
          features: this.features,
        });
        const pathGenerator = d3.geoPath<CellFeature>().projection(projection);
        return pathGenerator(d) || '';
      })
      .attr('stroke-width', (d: CellFeature) => {
        return d.properties.leiden === selectedCluster ? '3px' : '1px';
      })
      .attr('stroke', (d: CellFeature) => {
        return d.properties.leiden === selectedCluster ? '#000' : 'transparent';
      })
      // Remove mouseleave event to prevent resetting outline
      .on('mouseleave', null)
      .style('opacity', (d: CellFeature) => {
        return d.properties.leiden === selectedCluster ? 1.0 : 0.6;
      });
  }

  private getScaledPath(feature: CellFeature, scaleFactor: number): string {
    const coords = feature.geometry.coordinates[0];

    // Calculate centroid of the hexagon
    let centerX = 0,
      centerY = 0;
    coords.forEach((coord: number[]) => {
      centerX += coord[0];
      centerY += coord[1];
    });
    centerX /= coords.length;
    centerY /= coords.length;

    // Scale each coordinate outward from the center
    const scaledCoords = coords.map((coord: number[]) => {
      const dx = coord[0] - centerX;
      const dy = coord[1] - centerY;
      return [centerX + dx * scaleFactor, centerY + dy * scaleFactor];
    });

    // Create scaled geometry
    const scaledGeometry: CellGeometry = {
      type: 'Polygon',
      coordinates: [scaledCoords],
    };

    // Use path generator to convert to SVG path
    const projection = d3.geoIdentity().fitSize([1200, 1000], {
      type: 'FeatureCollection',
      features: this.features,
    });
    const pathGenerator = d3.geoPath().projection(projection);

    return pathGenerator(scaledGeometry) || '';
  }

  private resetClusterExtension(): void {
    const projection = d3.geoIdentity().fitSize([1200, 1000], {
      type: 'FeatureCollection',
      features: this.features,
    });
    const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

    this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => pathGenerator(d) || '')
      .attr('stroke-width', '1px')
      .attr('stroke', 'transparent')
      .style('opacity', 0.8);

    // Reinitialize the mouseleave event
    this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .on('mouseleave', (event, d) => this.mouseLeave(event, d));
  }

  private renderFootprintPlots(dataset: Dataset | null): void {
    this.footprintPlotUrls = (dataset?.footprint_list ?? []).map(relativePath =>
      this.sanitizer.bypassSecurityTrustResourceUrl(
        `${this.sessionService.apiUrl}/api/download/${relativePath}`
      )
    );
    // Load available motifs for the on-demand compute form.
    // Pass dataset_id so the backend can resolve the adata path even when
    // /read_adata has not been called yet (e.g. rescanned datasets).
    if (dataset) {
      const params = dataset.id ? `?dataset_id=${encodeURIComponent(dataset.id)}` : '';
      this.http.get<{ motifs: string[] }>(
        `${this.sessionService.apiUrl}/api/motifs${params}`,
        { withCredentials: true }
      ).subscribe({
        next: resp => { this.availableMotifs = resp.motifs ?? []; },
        error: () => { this.availableMotifs = []; }
      });
    } else {
      this.availableMotifs = [];
    }
  }

  public computeFootprint(): void {
    if (!this.footprintMotif) return;
    this.isComputingFootprint = true;
    this.footprintComputeError = '';
    const body = new FormData();
    body.append('motif', this.footprintMotif);
    body.append('cluster_by', this.footprintClusterBy);
    // Pass dataset_id for rescanned-dataset fallback (session may not have adata_path set)
    if (this.selectedDataset?.id) {
      body.append('dataset_id', this.selectedDataset.id);
    }
    this.http.post<{ footprint_url: string; relative_path: string }>(
      `${this.sessionService.apiUrl}/api/compute_footprint`,
      body,
      { withCredentials: true }
    ).subscribe({
      next: resp => {
        this.isComputingFootprint = false;
        const url = this.sanitizer.bypassSecurityTrustResourceUrl(
          `${this.sessionService.apiUrl}${resp.footprint_url}`
        );
        this.onDemandFootprintUrls = [...this.onDemandFootprintUrls, url];
      },
      error: err => {
        this.isComputingFootprint = false;
        this.footprintComputeError =
          err?.error?.detail ?? 'Footprint computation failed. Check the server logs.';
      }
    });
  }

    private renderNhoodHeatmap(): void {
      const leiden = this.selectedCell?.properties?.leiden;
      const clusterAnnotation = this.getLeidenClusterAnnotation(leiden);
      const enrichment = clusterAnnotation?.neighborhood_enrichment;

      if (!enrichment || !Array.isArray(enrichment)) return;

      const n = enrichment.length;
      const clusterLabels = Array.from({ length: n }, (_, i) => `Cluster ${i}`);

      const minValue = Math.min(...enrichment);
      const maxValue = Math.max(...enrichment);
      const normalized = (maxValue > minValue)
        ? enrichment.map((v: number) => (v - minValue) / (maxValue - minValue))
        : enrichment.map(() => 0);

      const data: Partial<Plotly.PlotData>[] = [
        {
          x: clusterLabels,
          y: normalized,
          type: 'bar',
          marker: { color: 'rgba(55, 128, 191, 0.7)' },
          name: `Cluster ${leiden} Neighborhood Enrichment`,
        }
      ];

      const layout = {
        margin: { t: 30, l: 60, r: 10, b: 40 },
        width: 300,
        height: 170,
        xaxis: {
          title: { text: 'Cluster' },
          automargin: true,
          tickfont: { size: 10 },
        },
        yaxis: {
          title: { text: 'Enrichment' },
          automargin: true,
          tickfont: { size: 10 },
        },
      };

      const container = document.getElementById('cluster-nhood-heatmap');
      if (!container) {
        console.error('Container cluster-nhood-heatmap not found for rendering heatmap');
        return;
      }

      Plotly.purge(container);
      Plotly.newPlot(container, data, layout, { displayModeBar: false});
    }

  public closeSidenav(): void {
    this.selectedCell = null;
    this.updateHexColors();
  }

  public updateCoOccurrenceTable(): void {
    if (this.selectedCluster === null) {
      this.coOccurrenceData = [];
      return;
  }

    const clusterAnnotation = this.getLeidenClusterAnnotation(this.selectedCluster);
    const coOccurrenceMatrix = clusterAnnotation?.co_occurrence;


    if (!Array.isArray(coOccurrenceMatrix)) {
        console.warn('No co-occurrence data found for cluster', this.selectedCluster);
        this.coOccurrenceData = [];
        return;
    }

    this.coOccurrenceData = [];
    console.log('Selected Cluster:', this.selectedCluster, 'Interval:', this.selectedInterval);
    console.log('Matrix structure:', coOccurrenceMatrix);
    try {
      for (let j = 0; j < this.clusterCount; j++) {
        if (
          Array.isArray(coOccurrenceMatrix[j]) &&
          typeof coOccurrenceMatrix[j][this.selectedInterval] === 'number'
        ) {
          this.coOccurrenceData.push(
            coOccurrenceMatrix[j][this.selectedInterval],
          );
        } else {
          this.coOccurrenceData.push(0);
        }
      }
    } catch (error) {
      console.error('Error extracting co-occurrence data:', error);
      this.coOccurrenceData = Array(this.clusterCount).fill(0);
    }

      // Calculate threshold for highlighting
      this.calculateCoOccurrenceThreshold();

      console.log(
          'Co-occurrence data for cluster',
          this.selectedCluster,
          'at interval',
          this.selectedInterval,
          ':',
          this.coOccurrenceData,
      );
    }

  private calculateCoOccurrenceThreshold(): void {
    const allValues = this.coOccurrenceData.flat().filter((val) => val > 0);
    if (allValues.length > 0) {
      allValues.sort((a, b) => a - b);
      const percentile75 = Math.floor(allValues.length * 0.75);
      this.coOccurrenceThreshold = allValues[percentile75] || 0.5;
    }
  }

  public getCoOccurrenceColor(value: number): string {
    if (value === 0) return '#f8f9fa';

    // Use the same color scale as hexagons for consistency
    // Exclude the same cluster (diagonal) from max calculation to better scale other values
    const coOccurrenceWithoutSameCluster = this.coOccurrenceData.filter(
      (_, index) => index !== this.selectedCluster,
    );
    const maxValue = Math.max(...coOccurrenceWithoutSameCluster);
    const intensity = Math.min(value / maxValue, 1);
    return this.continuousColorScale(intensity);
  }

  public getIntervalStats(): { min: number; max: number; avg: number } {
    const allValues = this.coOccurrenceData.flat().filter((val) => val > 0);
    if (allValues.length === 0) return { min: 0, max: 0, avg: 0 };

    const min = Math.min(...allValues);
    const max = Math.max(...allValues);
    const avg = allValues.reduce((sum, val) => sum + val, 0) / allValues.length;

    return { min, max, avg: Math.round(avg * 100) / 100 };
  }

  async getRegulatoryScoresforSpots(barcode: string) {
    this.sessionService.callWithSession(() =>
      this.http.get(
        `${this.sessionService.apiUrl}/obsm/regulatory_scores/cell/${barcode}`,
        { withCredentials: true },
      ),
    ).subscribe({
      next: (res) => {
        const rawData = res as { [scoreType: string]: { [element: string]: number } };

        const genie3Data: TableData = {};
        const spongeData: TableData = {};
        const genie3ElementsSet = new Set<string>();
        const spongeElementsSet = new Set<string>();

        // Separate the raw score groups by suffix (_genie3 or _sponge) for separate tables
        for (const [scoreType, scores] of Object.entries(rawData)) {
          if (scoreType.endsWith('_genie3')) {
            genie3Data[scoreType] = scores;
            Object.keys(scores).forEach(element => genie3ElementsSet.add(element));
          } else if (scoreType.endsWith('_sponge')) {
            spongeData[scoreType] = scores;
            Object.keys(scores).forEach(element => spongeElementsSet.add(element));
          }
        }
        this.genie3RawData = genie3Data;
        this.spongeRawData = spongeData;

        // Features/Columns assigned as the list of elements (TFs or Genes)
        this.genie3Elements = Array.from(genie3ElementsSet);
        this.spongeElements = Array.from(spongeElementsSet);

      },
      error: (err) => {
        this.genie3RawData = {};
        this.spongeRawData = {};
        this.genie3Elements = [];
        this.spongeElements = [];
        console.error(
          `[Backend] Failed to load regulatory scores for ${barcode}`,
          err,
        );
      }
    });
  }

  async fetchAndUpdate(columnName: string, index: string, updateCompare: boolean = false) {
    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/obsm/${columnName}/${index}`,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          const data = res as { [barcode: string]: any };

          if (this.features) {
            for (const feature of this.features) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                feature.properties[this.colorByProperty] = data[barcode];
              }
            }
          }
          console.log('Features:', this.features);
          if (updateCompare && this.compareMode && this.dataCompare?.features) {
            for (const feature of this.dataCompare.features) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                feature.properties[this.selectedCompareView] = data[barcode];
              }
            }
            console.log(`[Backend] Also updated compare view property '${this.selectedCompareView}' from obsm["${columnName}"][${index}]`);
          }
          console.log(`[Backend] Loaded adata.obsm["${columnName}][${index}]`);
          this.updateHexColors();
        },
        error: (err) =>
          console.error(
            `[Backend] Failed to load adata.obsm["${columnName}][${index}]`,
            err,
          ),
      });
  }

  public onRegulatoryScoreChange(): void {
    if (
      this.selectedRegulatoryScore?.endsWith('genie3') &&
      this.selectedGeneSetGenie3 &&
      this.selectedGeneSetSponge
    ) {
      this.fetchAndUpdate(
        this.selectedRegulatoryScore,
        this.selectedGeneSetGenie3,
      );
    } else if (
      this.selectedRegulatoryScore?.endsWith('sponge') &&
      this.selectedGeneSetSponge
    ) {
      this.fetchAndUpdate(
        this.selectedRegulatoryScore,
        this.selectedGeneSetSponge,
      );
    }
  }


  keyCompareByLabel = (a: KeyValue<string, unknown>, b: KeyValue<string, unknown>) => {
    return this.label(a.key).localeCompare(this.label(b.key), 'de', { sensitivity: 'base' });
  };

  private expandedProps = new Set<string>();

  isArray(v: any): v is any[] {
    return Array.isArray(v);
  }

  isNestedArray(v: any): v is any[][] {
    return Array.isArray(v) && v.length > 0 && v.every(row => Array.isArray(row) || this.looksLikeArrayString(row));
  }

  isNumberLike(v: unknown): v is number | string {
    return (typeof v === 'number' && Number.isFinite(v)) ||
      (typeof v === 'string' && v.trim() !== '' && Number.isFinite(+v));
  }

  toNumberLike(v: number | string): number {
    return typeof v === 'number' ? v : Number(v);
  }

  isNumericArray(arr: any): arr is (number | string)[] {
    return Array.isArray(arr) && arr.length > 0 && arr.every(x => this.isNumberLike(x));
  }

  isPrimitive(v: unknown): v is string | number | boolean | null {
    return v === null || ['string', 'number', 'boolean'].includes(typeof v as string);
  }

  getArrayStats(arr: (number | string)[]) {
    const nums = arr.map(x => Number(this.toNumberLike(x)));
    const min = Math.min(...nums);
    const max = Math.max(...nums);
    const avg = nums.reduce((a, b) => a + b, 0) / nums.length;
    return { min, max, avg };
  }

  formatValue(v: unknown): string {
    if (Array.isArray(v)) return v.join(', ');
    if (v && typeof v === 'object') return JSON.stringify(v, null, 2);
    return String(v);
  }

  label(key: string): string {
    return key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
  }

  toggleExpand(key: string) {
    if (this.expandedProps.has(key)) this.expandedProps.delete(key);
    else this.expandedProps.add(key);
  }

  isExpanded(key: string) {
    return this.expandedProps.has(key);
  }

  trackByIndex(index: number) { return index; }

  asArrayRow(row: any): any[] {
    if (Array.isArray(row)) return row;
    if (this.looksLikeArrayString(row)) {
      try {
        const parsed = JSON.parse(row as string);
        return Array.isArray(parsed) ? parsed : [row];
      } catch {
        return [row];
      }
    }
    return [row];
  }

  looksLikeArrayString(v: any): v is string {
    return typeof v === 'string' && /^\s*\[.*\]\s*$/.test(v);
  }

  toJsonCompact(obj: unknown, max = 120): string {
    try {
      const s = JSON.stringify(obj);
      return s.length > max ? s.slice(0, max) + '…' : s;
    } catch {
      return String(obj);
    }
  }


  shouldShowProperty(key: string): boolean {
    if (key == null) return true;
    const k = String(key).toLowerCase();
    return !this.hiddenPropKeys.has(k);
  }

  // ----- Dict/Object helpers -----
  isPlainObject(v: any): v is Record<string, any> {
    return v !== null && typeof v === 'object' && !Array.isArray(v);
  }

  objectKeyCount(obj: Record<string, any>): number {
    return Object.keys(obj).length;
  }

  objectEntries(obj: Record<string, any>): Array<{ key: string; value: any }> {
    return Object.keys(obj).sort().map(k => ({ key: k, value: obj[k] }));
  }

  prettyKey(k: string): string {
    return k.replace(/[_-]+/g, ' ')
      .replace(/\b\w/g, c => c.toUpperCase());
  }

  dictId(propLabel: string): string {
    return `DICT::${propLabel}`;
  }

  public selectViewTutorial(): void {
    const tour = new Shepherd.Tour({
      useModalOverlay: true,
      defaultStepOptions: {
        classes: 'shepherd-theme-custom'
      }
    });

    tour.addStep({
      id: 'info-btn-dropdown',
      attachTo: { element: '#color-by-property', on: 'left' },
      text: 'This button allows you to select the view you want to access. You can choose different regulatory scores, LIANA+ scores, Leiden Clustering and more. You can color by your obs columns in the view',
      buttons: [{ text: "Done", action: tour.complete }]
    });

    tour.start();
  }

  private renderLegend(containerName: string): void {
    const viewVariablesToUpdate = this.getViewVariablesToUpdate(containerName);

    viewVariablesToUpdate.svg.selectAll(`.${viewVariablesToUpdate.legendContainerName}`).remove();

    if (viewVariablesToUpdate.getLegendType() === 'continuous') {
      const [min, max] = viewVariablesToUpdate.getLegendDomain() as number[] || [0, 1];
      const legendX = 0;
      const legendY = 50;
      const width = 250;
      const height = 30;
      const fontSize = 24;
      const padding = 15;

      // Use standard <defs>
      const defs = viewVariablesToUpdate.svg.select('defs').empty()
        ? viewVariablesToUpdate.svg.append('defs')
        : viewVariablesToUpdate.svg.select('defs');

      defs.select(`#${viewVariablesToUpdate.legendGradientName}`).remove();

      const gradient = defs
        .append('linearGradient')
        .attr('id', viewVariablesToUpdate.legendGradientName)
        .attr('x1', '0%')
        .attr('x2', '100%')
        .attr('y1', '0%')
        .attr('y2', '0%');

      const numStops = 10;
      for (let i = 0; i <= numStops; i++) {
        const t = i / numStops;
        const value = min + t * (max - min);
        gradient
          .append('stop')
          .attr('offset', `${t * 100}%`)
          .attr('stop-color', viewVariablesToUpdate.continuous(value));
      }

      const legendG = viewVariablesToUpdate.svg
        .append('g')
        .attr('class', viewVariablesToUpdate.legendContainerName)
        .attr('transform', `translate(${legendX},${legendY})`);

      const titleText = this.translationService.translateSync(viewVariablesToUpdate.view);
      // fallback if translation returns empty
      const legendTitle = titleText && String(titleText).trim() ? titleText : this.label(this.selectedCompareView);

      // measure sizes using svg_compare
      const tempSvg = viewVariablesToUpdate.svg.append('g').style('opacity', 0);
      const titleWidth =
        tempSvg
          .append('text')
          .text(titleText)
          .style('font-size', `${fontSize}px`)
          .style('font-weight', 'bold')
          .node()
          ?.getBBox().width || 0;

      const minText = (min ?? 0).toFixed(2);
      const minWidth =
        tempSvg
          .append('text')
          .text(minText)
          .style('font-size', `${fontSize}px`)
          .node()
          ?.getBBox().width || 0;

      const maxText = (max ?? 0).toFixed(2);
      const maxWidth =
        tempSvg
          .append('text')
          .text(maxText)
          .style('font-size', `${fontSize}px`)
          .node()
          ?.getBBox().width || 0;

      tempSvg.remove();

      const textHeight = fontSize * 1.2;
      const requiredWidth = Math.max(width, titleWidth, minWidth + maxWidth + 20);
      const bgWidth = requiredWidth + padding * 2;
      const bgHeight = height + textHeight * 2 + padding * 3;

      legendG
        .append('rect')
        .attr('x', -padding)
        .attr('y', -padding - textHeight)
        .attr('width', bgWidth)
        .attr('height', bgHeight)
        .style('fill', 'rgba(255, 255, 255, 0.9)')
        .attr('stroke', '#ccc')
        .attr('stroke-width', 1)
        .attr('rx', 5);


      legendG
        .append('rect')
        .attr('x', (bgWidth - width) / 2 - padding)
        .attr('y', 0)
        .attr('width', width)
        .attr('height', height)
        .style('fill', `url(#${viewVariablesToUpdate.legendGradientName})`)
        .attr('stroke', '#ccc')
        .attr('stroke-width', 1)
        .attr('rx', 3);

      // Min label
      legendG
        .append('text')
        .attr('x', (bgWidth - width) / 2 - padding)
        .attr('y', height + textHeight)
        .attr('text-anchor', 'start')
        .style('font-size', `${fontSize}px`)
        .style('fill', '#333')
        .text(minText);

      // Max label
      legendG
        .append('text')
        .attr('x', (bgWidth - width) / 2 - padding + width)
        .attr('y', height + textHeight)
        .attr('text-anchor', 'end')
        .style('font-size', `${fontSize}px`)
        .style('fill', '#333')
        .text(maxText);

      // Title (compare) — position inside background with padding so it's not clipped
      const titleY = -padding + Math.round(fontSize / 2);
      legendG
        .append('text')
        .attr('x', bgWidth / 2 - padding)
        .attr('y', titleY)
        .attr('text-anchor', 'middle')
        .style('font-size', `${fontSize}px`)
        .style('font-weight', 'bold')
        .style('fill', '#333')
        .text(legendTitle);

    } else {
      // Categorical legend
      const categories = viewVariablesToUpdate.getLegendDomain() as string[] || [];
      categories.sort();
      const legendX = 0;
      const legendY = 50;
      const itemHeight = 40;
      const rectHeight = 20;
      const rectWidth = 30;
      const fontSize = 24;
      const titlePadding = 15;
      const padding = 15;

      // Measure using svg
      const tempSvg = viewVariablesToUpdate.svg.append('g').style('opacity', 0);
      const titleText = this.translationService.translateSync(viewVariablesToUpdate.view);
      const legendTitleCat = titleText && String(titleText).trim() ? titleText : this.label(viewVariablesToUpdate.view);
      const titleWidth =
        tempSvg
          .append('text')
          .text(legendTitleCat)
          .style('font-size', `${fontSize}px`)
          .style('font-weight', 'bold')
          .node()
          ?.getBBox().width || 0;

      const textNodes = tempSvg
        .selectAll('text')
        .data(categories)
        .enter()
        .append('text')
        .text((d) => d)
        .style('font-size', `${fontSize}px`);

      const maxTextWidth = categories.length
        ? Math.max(...textNodes.nodes().map((node) => (node as SVGGraphicsElement).getBBox().width))
        : 0;
      tempSvg.remove();

      const itemWidth = Math.max(200, maxTextWidth + rectWidth + 20, titleWidth + 40);
      const backgroundWidth = itemWidth + padding * 2;
      const categoryItemsHeight = categories.length * itemHeight;
      const titleHeight = fontSize * 1.2 + titlePadding;
      const backgroundHeight = titleHeight + categoryItemsHeight + padding * 2;

      const legendG = viewVariablesToUpdate.svg
        .append('g')
        .attr('class', viewVariablesToUpdate.legendContainerName)
        .attr('transform', `translate(${legendX},${legendY})`);

      // Background rectangle for categorical legend
      legendG
        .append('rect')
        .attr('x', -padding)
        .attr('y', -padding)
        .attr('width', backgroundWidth)
        .attr('height', backgroundHeight)
        .style('fill', 'rgba(255, 255, 255, 0.9)')
        .attr('stroke', '#ccc')
        .attr('stroke-width', 1)
        .attr('rx', 5);

      // Title (categorical)
      legendG
        .append('text')
        .attr('x', backgroundWidth / 2 - padding)
        .attr('y', titlePadding + fontSize / 2)
        .attr('dy', '0.35em')
        .attr('text-anchor', 'middle')
        .style('font-size', `${fontSize}px`)
        .style('font-weight', 'bold')
        .style('fill', '#333')
        .text(legendTitleCat);

      categories.forEach((cat, i) => {
        const yPosition = i * itemHeight + titleHeight;
        const legendItem = legendG.append('g').attr('transform', `translate(0, ${yPosition})`);
        const rectY = (itemHeight - rectHeight) / 2;
        legendItem
          .append('rect')
          .attr('y', rectY)
          .attr('width', rectWidth)
          .attr('height', rectHeight)
          .style('fill', viewVariablesToUpdate.ordinal(cat))
          .attr('stroke', '#333')
          .attr('stroke-width', 0.5)
          .attr('rx', 2);
        const textY = rectY + rectHeight / 2;
        legendItem
          .append('text')
          .attr('x', rectWidth + 10)
          .attr('y', textY)
          .attr('dy', '0.35em')
          .style('font-size', `${fontSize}px`)
          .style('fill', '#333')
          .text(cat);
        legendItem.append('title').text(cat);
      });
    }
  }

  private updateDetailAtScreenPos(screenX: number, screenY: number) {
    if (!this.isXenium || !this.currentPathGenerator) return;

    const half = this.detailSize / 2;

    const x0s = screenX - half;
    const x1s = screenX + half;
    const y0s = screenY - half;
    const y1s = screenY + half;

    const [x0d, y0d] = this.currentTransform.invert([x0s, y0s]);
    const [x1d, y1d] = this.currentTransform.invert([x1s, y1s]);

    const centerX = (x0d + x1d) / 2;
    const centerY = (y0d + y1d) / 2;

    const [targetX, targetY] = this.currentTransform.invert([screenX, screenY]);

    const localScale = 6;

    this.detailLayer.attr(
      'transform',
      `translate(${targetX},${targetY}) scale(${localScale}) translate(${-centerX},${-centerY})`
    );

    const subset = this.fullFeatures.filter((f) => {
      const c = (f.properties as any).__centroidProjected;
      if (!c || c.length < 2) return false;

      const [x, y] = c as [number, number];
      return x >= x0d && x <= x1d && y >= y0d && y <= y1d;
    });

    this.detailLayer
        .style('cursor', 'pointer')
        .selectAll<SVGPathElement, CellFeature>('path')
        .data(subset)
        .join('path')
        .attr('d', (d: CellFeature) => this.currentPathGenerator!(d) || '')
        .attr('fill', (d: CellFeature) => {
            const value = this.leidenCentralityProps.includes(this.colorByProperty)
              ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[this.colorByProperty]
              : d.properties?.[this.colorByProperty];
            if (this.currentLegendType === 'categorical') {
                return this.colorScale(String(value));
            } else {
                const num = this.toNumber(value);
                return Number.isFinite(num)
                    ? this.continuousColorScale(num)
                    : '#ccc';
            }
        })
        .attr('stroke', '#fff')
        .attr('stroke-width', 0.4)
        .style('opacity', 1)
        .on('click', (event, d) => this.openSidenav(event, d));
}


  private initDetailWindow() {
    // alte ClipPaths/Frames entfernen (falls vorhanden)
    this.svg.select("#detail-frame").remove();

    // Clip-Fenster (immer in <defs>)
    const defs = this.svg.select('defs').empty()
      ? this.svg.append('defs')
      : this.svg.select('defs');

    defs.select("#detail-clip").remove();

    defs.append("clipPath")
      .attr("id", "detail-clip")
      .append("rect")
      .attr("id", "detail-window")
      .attr("x", 30)
      .attr("y", 30)
      .attr("width", this.detailSize)
      .attr("height", this.detailSize);


    // Weißer Hintergrund nur im Detail-Layer (Xenium)
    this.detailLayer
      .selectAll(".detail-bg")
      .remove();

    this.detailLayer
      .append("rect")
      .attr("class", "detail-bg")
      .attr("x", -10000)
      .attr("y", -10000)
      .attr("width", 20000)
      .attr("height", 20000)
      .attr("fill", "#ffffff")
      .attr("pointer-events", "none");

    // Rahmen des Fensters
    this.svg.append("rect")
      .attr("id", "detail-frame")
      .attr("x", 30)
      .attr("y", 30)
      .attr("width", this.detailSize)
      .attr("height", this.detailSize)
      .attr("fill", "none")
      .attr("stroke", "#111")
      .attr("stroke-width", 1)
      .style("pointer-events", "none");
  }

  private bindDetailWindowInteractions(): void {
    if (!this.svg) return;

    // Alte Handler entfernen, damit bei Reloads nichts doppelt gebunden wird
    this.svg.on('click.detail', null);
    this.svg.on('dblclick.detail', null);
    this.svg.on('contextmenu.detail', null);

    // Doppelklick = explizit repositionieren
    this.svg.on('dblclick.detail', (event: MouseEvent) => {
      event.preventDefault();
      event.stopPropagation();

      const [mx, my] = d3.pointer(event, this.svg.node());
      this.setDetailWindow(mx, my);
    });

    // Rechtsklick = entfernen
    this.svg.on('contextmenu.detail', (event: MouseEvent) => {
      event.preventDefault();
      event.stopPropagation();
      this.clearDetailWindow();
    });

    // ESC = schließen
    if (this.keydownHandler) {
      window.removeEventListener('keydown', this.keydownHandler);
    }

    this.keydownHandler = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        this.clearDetailWindow();
      }
    };

    window.addEventListener('keydown', this.keydownHandler);
  }

  private setDetailWindow(screenX: number, screenY: number): void {
    this.detailVisible = true;
    this.detailScreenPos = { x: screenX, y: screenY };

    this.showDetailWindowAt(screenX, screenY);
    this.updateDetailAtScreenPos(screenX, screenY);
  }

  private clearDetailWindow(): void {
    this.detailVisible = false;
    this.detailScreenPos = null;
    this.hideDetailWindow();

    // Optional: Detail-Layer Inhalte leeren
    this.detailLayer.selectAll('path').remove();
  }

  private showDetailWindowAt(screenX: number, screenY: number): void {
    this.svg.select('#detail-window')
      .attr('x', screenX - this.detailSize / 2)
      .attr('y', screenY - this.detailSize / 2)
      .style('display', null);

    this.svg.select('#detail-frame')
      .attr('x', screenX - this.detailSize / 2)
      .attr('y', screenY - this.detailSize / 2)
      .style('display', null);

    this.detailLayer.style('display', null);
  }

  private hideDetailWindow(): void {
    this.svg.select('#detail-window').style('display', 'none');
    this.svg.select('#detail-frame').style('display', 'none');
    this.detailLayer.style('display', 'none');
  }
}




interface CellGeometry {
  type: 'Polygon';
  coordinates: number[][][];
}

interface CellProperties {
    barcode: string;
    centroid: [number, number] | [];
    cell_type: string;
    leiden: number;
    color: string;
    aucell_genie3: { [key: string]: number };
    aucell_sponge: { [key: string]: number };
    [key: string]: string | number | number[] | [] | undefined | { [key: string]: any };
}

interface CellFeature {
  type: 'Feature';
  geometry: CellGeometry;
  properties: CellProperties;
}

interface GeoJsonData {
  type: 'FeatureCollection';
  features: CellFeature[];
  meta?: { [key: string]: any };
}

interface genie3RegGraphConnection {
  source: string;
  target: string;
  weight: number;
}

interface spongeRegGraphConnection {
  source: string;
  target: string;
  p_adjusted: number;
}

interface TableData {
  [columnHeader: string]: { [rowHeader: string]: string | number };
}
