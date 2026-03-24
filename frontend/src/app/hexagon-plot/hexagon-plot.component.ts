import { Component, OnInit, OnDestroy, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { KeyValue } from '@angular/common';
import { Data, Router } from '@angular/router';
import { firstValueFrom, map, Observable, Subject, Subscription, takeUntil } from 'rxjs';

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
import { MatSelect, MatSelectTrigger } from '@angular/material/select';
import { MatExpansionModule } from '@angular/material/expansion';
import { MatTableModule } from '@angular/material/table';
import { MatDividerModule } from '@angular/material/divider';
import { MatTabsModule } from '@angular/material/tabs';
import { MatTabChangeEvent, MatTabGroup } from '@angular/material/tabs';
import { DomSanitizer, SafeResourceUrl } from '@angular/platform-browser';
import { MatInputModule } from '@angular/material/input';
import { MatCheckboxModule } from '@angular/material/checkbox';



@Component({
  selector: 'app-hexagon-plot',
  imports: [CommonModule, FormsModule, FilterableTableComponent, TranslatePipe, MatButtonModule, MatIconModule, MatTooltipModule, MatDialogModule, MatProgressSpinnerModule, MatOptgroup, MatFormField, MatLabel, MatOption, MatSelect, MatSelectTrigger, MatExpansionModule, MatTableModule, MatDividerModule, MatTabsModule, MatInputModule, MatCheckboxModule],
  standalone: true,
  templateUrl: './hexagon-plot.component.html',
  styleUrls: ['./hexagon-plot.component.scss'],
})
export class HexagonPlotComponent implements OnInit, OnDestroy, AfterViewInit {
  @ViewChild('aucell_graph_genie3', { static: false }) aucellGraphGenie3Element?: ElementRef<HTMLElement>;
  @ViewChild('aucell_graph_sponge', { static: false }) aucellGraphSpongeElement?: ElementRef<HTMLElement>;
  @ViewChild(MatTabGroup, { static: false }) tabGroup?: MatTabGroup;
  @ViewChild('dgeaHeatmap', { static: false }) dgeaHeatmapElement!: ElementRef<HTMLElement>;
  private _resizeHandler: any = null;
  // Removed manual subscription tracker - all subscriptions use takeUntil now
  footprintPlotUrls: SafeResourceUrl[] = [];
  onDemandFootprintUrls: SafeResourceUrl[] = [];
  availableMotifs: string[] = [];
  availableCellTypes: string[] = [];
  footprintMotifs: string[] = [];
  motifSearchQuery: string = '';
  footprintClusterBy: string = 'cell_type';
  isComputingFootprint: boolean = false;
  footprintComputeError: string = '';

  get filteredMotifs(): string[] {
    const q = this.motifSearchQuery.trim().toLowerCase();
    return q ? this.availableMotifs.filter(m => m.toLowerCase().includes(q)) : this.availableMotifs;
  }

  builtinDatasets$: Observable<Dataset[]>;
  uploadedDatasets$: Observable<Dataset[]>;
  tangram_Datasets$: Observable<Dataset[]>;
  selectedDataset: Dataset | null = null;
  selectedDatasetCompare: Dataset | null = null;

  compareMode: boolean = false;

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
  public compareDataPath = DEFAULT_PATHS.hexagonPath;
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


  public selectedCell: CellFeature | null = null;
  public selectedCellCompare: CellFeature | null = null;
  public selectedCluster: number | null = null;
  public selectedClusterCompare: number | null = null;

  public selectedView = 'regulatory_scores';
  public selectedCompareView: string = 'regulatory_scores';


  public selectedGeneSetGenie3: string | null = null;
  public selectedGeneSetGenie3Compare: string | null = null;
  public selectedGeneSetSpongeCompare: string | null = null;
  public selectedGeneSetSponge: string | null = null;
  public selectedRegulatoryScore: string | null = null;
  public selectedRegulatoryScoreCompare: string | null = null;
  public selectedGeneExpressionMain: string | null = null;
  public selectedGeneExpressionCompare: string | null = null;
  private regulatoryObsmKeysMain: string[] = [];
  private regulatoryObsmKeysCompare: string[] = [];

  // Data sources for the two tables
  public genie3RawData: TableData = {};
  public spongeRawData: TableData = {};
  public genie3RawDataCompare: TableData = {};
  public spongeRawDataCompare: TableData = {};

  // Column lists for the two tables
  public genie3Elements: string[] = [];
  public spongeElements: string[] = [];
  public genie3ElementsCompare: string[] = [];
  public spongeElementsCompare: string[] = [];

  private previousGeneSetGenie3: string | null = null;
  private previousGeneSetSponge: string | null = null;
  private previousGeneSetGenie3Compare: string | null = null;
  private previousGeneSetSpongeCompare: string | null = null;
  private requestTokens: { [key: string]: number } = {};
  private regulatoryFetchTokenMain: number = 0;
  private regulatoryFetchTokenCompare: number = 0;
  private activeRegulatoryTabSource: 'main' | 'compare' = 'main';
  private geneDomainToken: number = 0;
  private geneDomainCache = new Map<string, { min: number; max: number; expiresAt: number }>();
  private geneDomainCacheTtlMs: number = 10 * 60 * 1000;
  private sharedGeneExpressionDomain: { min: number; max: number } | null = null;
  private sharedGeneExpressionContextKey: string | null = null;
  public dgeaReady: boolean = false;

  public selectedInterval: number = 0;
  public selectedIntervalCompare: number = 0;

  public features: CellFeature[] = [];
  public compareFeatures: CellFeature[] = [];

  public meta: { [key: string]: any } = {};
  public metaCompare: { [key: string]: any } = {};

  // Selected groups for the DGEA comparison (bound to the dropdowns)
  public selectedDgeaObsCol: string = 'cell_type';
  public selectedDgeaGroup1: string | null = null;
  public selectedDgeaGroup2: string | null = null;
  private hiddenPropKeys = new Set<string>([]);
  public dgeaVsAll: boolean = false;


  public clusterCells: CellFeature[] = [];
  public clusterCellsCompare: CellFeature[] = [];
  public clusterCellTypes: {
    type: string;
    count: number;
    percentage: string;
  }[] = [];
  public clusterCellTypesCompare: {
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
  public clusterCentralityAvgCompare: {
    degree_centrality: number;
    average_clustering: number;
    closeness_centrality: number;
  } = {
      degree_centrality: 0,
      average_clustering: 0,
      closeness_centrality: 0,
    };


  public compareClusterCells: CellFeature[] = [];
  public compareClusterCellTypes: {
    type: string;
    count: number;
    percentage: string;
  }[] = [];
  public compareClusterCentralityAvg: {
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

  public genie3NetworkCompare: genie3RegGraphConnection[] = [];
  public spongeNetworkCompare: spongeRegGraphConnection[] = [];

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

  public genie3WeightCutoffCompare: number = 0.5;
  public spongePValueCutoffCompare: number = 0.05;
  public genie3MinEdgesCompare: number = 25;
  public spongeMinEdgesCompare: number = 25;
  public spongeSliderDataCompare: { step: number; min_border: number; max_border: number; default_value: number } | null = null;
  public genie3SliderDataCompare: { step: number; min_border: number; max_border: number; default_value: number } | null = null;


  // Loading screen trackers
  public isLoadingHexagons: boolean = true;
  public isLoadingCompare: boolean = false;
  public isLoadingSponge: boolean = false;
  public isLoadingGenie3: boolean = false;
  public isLoadingGenie3Compare: boolean = false;
  public isLoadingSpongeCompare: boolean = false;

  public isLoadingRegulatoryScores: boolean = false; // Track regulatory scores fetches during init
  public isLoadingRegulatoryScoresCompare: boolean = false; // Track regulatory scores fetches during init for compare view
  public isAppInitializing: boolean = true; // Only for initial app load

  // Returns true if any initial data loading is in progress (ONLY for first time)
  public get isInitializing(): boolean {
    return this.isAppInitializing && (this.isLoadingHexagons || this.isLoadingGenie3 || this.isLoadingSponge || this.isLoadingRegulatoryScores);
  }

  // Co-occurrence table
  public coOccurrenceData: number[] = [];
  public coOccurrenceDataCompare: number[] = [];
  public coOccurrenceColumns: string[] = [];
  public coOccurrenceThreshold: number = 0.5;
  public coOccurrenceThresholdCompare: number = 0.5;
  public maxInterval: number = 49;
  public clusterCount: number = 10;

  // Compare Co-occurrence table
  public compareCoOccurrenceData: number[] = [];
  public compareCoOccurrenceColumns: string[] = [];
  public compareCoOccurrenceThreshold: number = 0.5;
  public maxIntervalCompare: number = 49;
  public clusterCountCompare: number = 10;

  public colorableProperties: string[] = [
    'cell_type',
    'leiden',
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];

  public colorablePropertiesCompare: string[] = [
    'cell_type',
    'leiden',
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];

  public geneSetsGenie3Compare: { [regulator: string]: string[] } = {};
  public geneSetsSpongeCompare: { [regulator: string]: string[] } = {};

  // keep track of which properties actually have any data available; used to
  // disable/grey‑out menu items and tabs instead of hiding them abruptly.
  public propertyAvailability: { [prop: string]: boolean } = {};
  public propertyAvailabilityCompare: { [prop: string]: boolean } = {};

  public leidenCentralityProps = [
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];

  public groupedProperties: { key: string; value: string[] }[] | null = null;
  public groupedPropertiesCompare: { key: string; value: string[] }[] | null = null;

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
  private continuousColorScale = d3.scaleSequential(d3.interpolateYlOrRd).clamp(true);
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
  private continuousColorScaleCompare = d3.scaleSequential(d3.interpolateYlOrRd).clamp(true);
  public currentLegendDomain: any[] = [];
  public currentLegendType: 'continuous' | 'categorical' = 'categorical';

  // Comparison

  public currentCompareLegendType: 'continuous' | 'categorical' = 'categorical';
  public currentLegendDomainCompare: any[] = [];


  public repaintBothViews(): void {
    this.updateHexColors();
    if (this.compareMode) {
      this.updateHexColors('#hexbin-compare');
    }
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

    // --- MAIN DATASET LISTENER ---
    this.datasetService.selectedDataset$
      .pipe(takeUntil(this.destroy$))
      .subscribe(dataset => {
        this.selectedDataset = dataset;
        if (dataset) {
          this.dataPath = dataset.geojson_path || '';
          this.features = [];
          this.meta = {};
          this.updatePathsFromDataset(dataset, false);
          // Note: updatePathsFromDataset triggers the pathsService listener,
          // which handles the actual D3 loading for the main view.
        } else {
          this.pathsService.updatePaths({ adataPath: undefined }, false);
          this.dataPath = '';
        }
      });

    // --- COMPARE DATASET LISTENER ---
    this.datasetService.selectedDatasetCompare$
      .pipe(takeUntil(this.destroy$))
      .subscribe(dataset => {
        this.selectedDatasetCompare = dataset;
        if (dataset && this.compareMode) {
          // Prepare the state
          this.compareDataPath = dataset.geojson_path || '';
          this.compareFeatures = [];
          this.metaCompare = {};
          this.regulatoryObsmKeysCompare = [];

          this.updatePathsFromDataset(dataset, true);
          this.refreshCompareRegulatoryAvailability();

          // Trigger the reload
          this.reloadComparisonView();
        } else if (!dataset && this.compareMode) {
          this.pathsService.updatePaths({ adataPath: undefined }, true);
          this.compareDataPath = '';
          d3.select('#hexbin-compare').selectAll('svg').remove();
        }
      });

    this.isLoadingHexagons = true;
    // Subscribe to path changes with proper cleanup
    this.pathsService.paths$
      .pipe(takeUntil(this.destroy$))
      .subscribe(paths => {
        if (!paths) return;

        const hexagonPath = paths.hexagonPath || DEFAULT_PATHS.hexagonPath;

        if (hexagonPath) {
          this.isLoadingHexagons = true; // 1. Turn on spinner
          this.dataPath = hexagonPath;

          d3.select('#hexbin').selectAll('svg').remove();
          this.footprintPlotUrls = [];
          this.onDemandFootprintUrls = [];
          this.availableMotifs = [];
          this.availableCellTypes = [];
          this.footprintMotifs = [];
          this.motifSearchQuery = '';
          this.footprintComputeError = '';

          // 2. Yield to browser to paint spinner, then do the heavy lifting
          setTimeout(() => {
            this.createHexagonPlot();
            this.loadAndRenderData(this.dataPath);
            this.renderFootprintPlots(this.selectedDataset);
          }, 50);
        } else {
          console.warn('✗ No hexagon path available');
        }
      });
  }

  ngOnDestroy(): void {

    this.destroy$.next();
    this.destroy$.complete();

    // Remove resize listener if added
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
    this._resizeHandler = () => this.updateTabPagination();
    window.addEventListener('resize', this._resizeHandler);
  }

  onDatasetSelected(dataset: Dataset | null, isCompare: boolean = false): void {
    if (isCompare) {
      // Just emit to the service. The listener in ngOnInit will catch it.
      this.datasetService.selectDatasetCompare(dataset);
    } else {
      // Just emit to the service.
      this.datasetService.selectDataset(dataset);
    }
  }

  private updateGraphWidths(): void {
    const genie3El = this.aucellGraphGenie3Element?.nativeElement;
    const spongeEl = this.aucellGraphSpongeElement?.nativeElement;

    if (genie3El) {
      this.genie3Width = genie3El.clientWidth || 600;
    }
    if (spongeEl) {
      this.spongeWidth = spongeEl.clientWidth || 600;
    }
  }

  public onRegulatoryGraphsPanelOpened(): void {
    // Wait for lazy expansion-panel content to be attached to the DOM.
    setTimeout(() => {
      this.updateGraphWidths();

      if (this.selectedGeneSetGenie3) {
        this.updateSubgraphGenie3();
      }

      if (this.selectedGeneSetSponge) {
        this.updateSubgraphSponge();
      }
    }, 0);
  }

  // Handle tangram dataset selection - use tangram_adata_path if available
  onDatasetTangramSelected(dataset: Dataset | null, isCompare: boolean = false): void {
    if (!dataset || !dataset.tangram_adata_path) return;

    // Create a modified dataset with tangram path
    const tangramDataset: Dataset = {
      ...dataset,
      adata_path: dataset.tangram_adata_path
    };

    if (isCompare) {
      this.selectedDatasetCompare = tangramDataset;
      this.datasetService.selectDatasetCompare(tangramDataset);
      this.updatePathsFromDataset(tangramDataset, true);
      this.reloadComparisonView();
    } else {
      this.selectedDataset = tangramDataset;
      this.datasetService.selectDataset(tangramDataset);
      this.updatePathsFromDataset(tangramDataset, false);
      this.reloadHexagons();
    }

    this.refreshSharedGeneExpressionDomain();
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

      // Yield to the browser's paint cycle so the spinner appears
      setTimeout(() => {
        this.loadAndRenderData(this.compareDataPath, true);
      }, 50);
    }
  }

  private updatePathsFromDataset(dataset: Dataset, isCompare: boolean = false): void {
    this.pathsService.updatePaths({
      adataPath: dataset.adata_path,
      genieFiltPath: dataset.genie_network_path,
      spongeFiltPath: dataset.sponge_network_path,
      hexagonPath: dataset.geojson_path,
    }, isCompare);
  }

  // allow nested tables like for differential motif activity view
  public asTableData(value: unknown): {
    [col: string]: { [index: string]: string | number }
  } | string[] {
    return value as { [col: string]: { [index: string]: string | number } } | string[];
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

  public onTabColorChange(newView: string, compare: boolean = false): void {
    if (compare) {
      if (this.selectedCompareView !== newView) {
        this.selectedCompareView = newView;
        this.onColorbyPropertyChange(true);
      } else {
        this.updateHexColors('#hexbin-compare');
      }
    } else {
      if (this.selectedView !== newView) {
        this.selectedView = newView;
        this.onColorbyPropertyChange(false);
      } else {
        this.updateHexColors();
      }
    }
  }

  public onTabChange(event: MatTabChangeEvent, compare: boolean = false): void {
    let newView: string | null = null;
    const tabLabel = event.tab.textLabel;

    // Handle info tabs that need async rendering (main view only)
    if (!compare && (tabLabel === 'Cluster Information' || tabLabel === 'Cell Information')) {
      setTimeout(() => this.renderNhoodHeatmap(), 300);
      return;
    }

    if (!compare && tabLabel === 'DGEA') {
      this.dgeaReady = !!this.meta?.['dgea']?.[this.selectedDgeaObsCol];
      if (this.dgeaReady) {
        this.initDgeaSelection();
        setTimeout(() => this.renderDgeaHeatmap(), 100);
      }
      return;
    }

    // Map tab labels to view keys
    const tabMap: { [label: string]: string } = compare
      ? {
        'Compare - Regulatory Scores': 'regulatory_scores',
        'Compare - Co-occurence': 'leiden',
        'Compare - Gene Expression': 'gene_expression',
        'Compare - Ligand-Receptor Relationships': 'ligand_receptor_relationships',
        'Compare - Cell Composition TF Activity': 'cell_comp_tf_activity_similarity',
        'Compare - TF Activity': 'tf_activity',
        'Compare - Pathway Activity': 'pathway_activity',
      }
      : {
        'Regulatory Scores': 'regulatory_scores',
        'Co-occurence': 'leiden',
        'Gene Expression': 'gene_expression',
        'TF Activity': 'tf_activity',
        'Pathway Activity': 'pathway_activity',
        'Cell Composition TF Activity': 'cell_comp_tf_activity_similarity',
        'Ligand-Receptor Relationships': 'ligand_receptor_relationships',
        "ChromVar spatial correlation : Moran's I / Geary's C": 'chromvar_total_sum',
        'Differential Motif Activity': 'cell_type',
        'Footprints': 'cell_type',
      };

    newView = tabMap[tabLabel] || null;

    console.log('[Tab Change]', compare ? '[Compare]' : '[Main]', 'newView=', newView);

    if (newView) {
      this.onTabColorChange(newView, compare);
    }
  }

  /**
   * Check if regulatory scores data is available
   */
  private hasRegulatoryScoresData(scope: 'main' | 'compare' = 'main'): boolean {
    const meta = scope === 'compare' ? this.metaCompare : this.meta;
    if (!meta) return false;

    const hasGenie3 = meta['global_regulatory_scores_genie3'] &&
      Object.keys(meta['global_regulatory_scores_genie3']).length > 0;
    const hasSponge = meta['global_regulatory_scores_sponge'] &&
      Object.keys(meta['global_regulatory_scores_sponge']).length > 0;

    const hasRegulatoryObsm = (
      scope === 'compare' ? this.regulatoryObsmKeysCompare : this.regulatoryObsmKeysMain
    ).length > 0;

    return (hasGenie3 || hasSponge) && hasRegulatoryObsm;
  }

  private hasGeneExpressionData(scope: 'main' | 'compare' = 'main'): boolean {
    const meta = scope === 'compare' ? this.metaCompare : this.meta;
    if (!meta) return false;

    const hasMoranI = meta['moranI'] && Object.keys(meta['moranI']).length > 0;
    const hasGearyC = meta['gearyC'] && Object.keys(meta['gearyC']).length > 0;

    return !!(hasMoranI || hasGearyC);
  }

  private async refreshCompareRegulatoryAvailability(): Promise<void> {
    const compareDatasetId = this.selectedDatasetCompare?.id || this.selectedDataset?.id;
    this.regulatoryObsmKeysCompare = await this.fetchRegulatoryObsmKeys(compareDatasetId);
  }

  private async fetchRegulatoryObsmKeys(datasetId?: string): Promise<string[]> {
    const query = datasetId ? `?dataset_id=${encodeURIComponent(datasetId)}` : '';
    try {
      const res = await firstValueFrom(
        this.sessionService.callWithSession(() =>
          this.http.get<{ regulatory_obsm_keys?: string[] }>(
            `${this.sessionService.apiUrl}/api/obsm_tables${query}`,
            { withCredentials: true },
          )
        )
      );
      return Array.isArray(res?.regulatory_obsm_keys) ? res.regulatory_obsm_keys : [];
    } catch (e) {
      console.warn('[Regulatory availability] Could not fetch obsm keys', e);
      return [];
    }
  }

  /**
   * Returns true if any data exists for the given property/column.  This is used
   * by the template to disable tabs and dropdown options when information
   * hasn't been computed or loaded yet.
   */
  public propertyAvailable(prop: string, scope: 'main' | 'compare' = 'main'): boolean {
    if (prop === 'gene_expression') return this.hasGeneExpressionData(scope)
    if (prop === 'regulatory_scores') return this.hasRegulatoryScoresData(scope);
    if (prop === 'co_occurrence') return this.hasCoOccurrenceData(scope);


    const map = scope === 'compare' ? this.propertyAvailabilityCompare : this.propertyAvailability;

    // Use the carefully calculated map if it exists for this property
    if (map && prop in map) {
      return map[prop];
    }

    return false;
  }

  private hasCoOccurrenceData(scope: 'main' | 'compare' = 'main'): boolean {
    const meta = scope === 'compare' ? this.metaCompare : this.meta;

    if (!meta || !meta['leiden_cluster_annotations']) return false;

    const annotations = meta['leiden_cluster_annotations'];

    // Check if at least one cluster has a valid 2D co_occurrence array
    return Object.values(annotations).some((ann: any) =>
      ann &&
      Array.isArray(ann.co_occurrence) &&
      ann.co_occurrence.length > 0 &&
      Array.isArray(ann.co_occurrence[0]) // verify it's a 2D array
    );
  }

  private updatePropertyAvailability(scope: 'main' | 'compare' = 'main') {
    const features = scope === 'main' ? this.features : this.compareFeatures;
    const map: { [prop: string]: boolean } = {};
    const currentColorableProps = scope === 'main' ? this.colorableProperties : this.colorablePropertiesCompare;
    if (features && features.length > 0) {
      const keys = Object.keys(features[0].properties || {});
      for (const prop of currentColorableProps) {
        if (prop === 'regulatory_scores') {
          map[prop] = this.hasRegulatoryScoresData(scope);
        } else {
          map[prop] = keys.includes(prop);
        }
      }
      // Add all other keys as available
      for (const key of keys) {
        map[key] = true;
      }
    }
    if (scope === 'main') {
      this.propertyAvailability = map;
    } else {
      this.propertyAvailabilityCompare = map;
    }
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

  public onCompareMode(): void {
    this.compareMode = !this.compareMode;

    if (this.compareMode) {
      this.syncCompareSelectionsFromMain();
      this.isLoadingCompare = true;
      this.regulatoryObsmKeysCompare = [...this.regulatoryObsmKeysMain];

      // Emit the main dataset to the compare stream.
      // This will automatically trigger the listener in ngOnInit!
      if (this.selectedDataset) {
        this.datasetService.selectDatasetCompare(this.selectedDataset);
      }
      this.refreshSharedGeneExpressionDomain();

    } else {
      // Teardown logic stays exactly the same
      d3.select('#hexbin-compare').selectAll('*').remove();
      try { d3.select('#hexbin-compare').selectAll('svg').remove(); } catch { }
      this.svg_compare = null as any;
      this.g_compare = null as any;

      try {
        this.currentLegendDomainCompare = [];
        this.currentCompareLegendType = 'categorical';
        if (this.colorScaleCompare && typeof this.colorScaleCompare.domain === 'function') {
          this.colorScaleCompare.domain([] as any);
        }
        if (this.continuousColorScaleCompare && typeof this.continuousColorScaleCompare.domain === 'function') {
          this.continuousColorScaleCompare.domain([0, 1]);
        }
      } catch (e) { }

      // Clear the service state
      this.datasetService.selectDatasetCompare(null);
      this.refreshSharedGeneExpressionDomain();
    }
  }

  private syncCompareSelectionsFromMain(): void {
    this.selectedCompareView = this.selectedView;

    if (this.selectedView === 'gene_expression') {
      this.selectedGeneExpressionCompare = this.selectedGeneExpressionMain;
    }

    if (this.selectedView === 'regulatory_scores') {
      if (
        this.selectedRegulatoryScore?.endsWith('genie3') &&
        this.selectedGeneSetGenie3
      ) {
        this.fetchAndUpdate(this.selectedRegulatoryScore, this.selectedGeneSetGenie3, {
          updateMain: true,
          updateCompare: true,
        });
      } else if (
        this.selectedRegulatoryScore?.endsWith('sponge') &&
        this.selectedGeneSetSponge
      ) {
        this.fetchAndUpdate(this.selectedRegulatoryScore, this.selectedGeneSetSponge, {
          updateMain: true,
          updateCompare: true,
        });
      }
    }
  }

  // Selected groups for the DGEA comparison (bound to the dropdowns)
  getSelectedDgeaHeatmap(): any | null {
    const cmp = this.getSelectedDgeaComparison();
    if (!cmp || cmp['skipped']) return null;
    return cmp['heatmap_context'] ?? null;
  }

  // Returns all available cell type levels for the dropdown selectors
  getDgeaLevels(): string[] {
    return this.meta?.['dgea']?.[this.selectedDgeaObsCol]?.['levels'] ?? [];
  }

  // Returns the map of all DGEA comparisons
  getDgeaComparisonMap(): { [key: string]: any } {
    return this.meta?.['dgea']?.[this.selectedDgeaObsCol]?.['comparisons'] ?? {};
  }

  // Re-render the DGEA heatmap when the user changes the group selections
  public onDgeaSelectionChange(): void {
    setTimeout(() => this.renderDgeaHeatmap(), 0);
  }

  getDgeaObsCols(): string[] {
    return Object.keys(this.meta?.['dgea'] ?? {});
  }

  public onDgeaObsColChange(): void {
    this.selectedDgeaGroup1 = null;
    this.selectedDgeaGroup2 = null;
    this.dgeaVsAll = false;
    this.dgeaReady = !!this.meta?.['dgea']?.[this.selectedDgeaObsCol];
    this.initDgeaSelection();
    setTimeout(() => this.renderDgeaHeatmap(), 0);
  }

  getDgeaObsColLabel(col: string): string {
    if (col === 'cell_type') return 'Cell type';
    if (col === 'leiden') return 'Leiden';
    return col;
  }

  public hasDgeaData(): boolean {
    return !!this.meta?.['dgea'] && Object.keys(this.meta['dgea']).length > 0;
  }


  private getLeidenClusterAnnotation(clusterId: number | null | undefined, compare: boolean = false): any | null {
    if (clusterId === null || clusterId === undefined) return null;
    const meta = compare ? this.metaCompare : this.meta;
    return meta?.['leiden_cluster_annotations']?.[String(clusterId)] ?? null;
  }

  // Create comparison ID matching the backend JSON format
  private makeComparisonId(group1: string, group2: string): string {
    const safe = (x: string) => x.replace(/[^A-Za-z0-9]+/g, '_');
    return `${safe(group1)}__vs__${safe(group2)}`;
  }

  public shownGeneOnPlot: string | null = null;

  public onDgeaGeneSelected(event: { gene: string; action: string }): void {
    if (event.action === 'show_on_plot') {
      this.showDgeaGeneOnMainPlot(event.gene);
    }
  }

  public showDgeaGeneOnMainPlot(gene: string): void {
    this.shownGeneOnPlot = gene;
    this.selectedView = 'gene_expression';
    this.fetchAndUpdate('gene_expression', gene);
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

    const x = (hm.groups as any[]).map(g => String(g));
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

    const xAxisTitle = this.getDgeaObsColLabel(this.selectedDgeaObsCol);

    const layout: Partial<Plotly.Layout> = {
      margin: { t: 30, l: 140, r: 20, b: 100 },
      height: Math.max(420, y.length * 22),
      xaxis: {
        title: { text: xAxisTitle },
        tickangle: -45,
        automargin: true,
        type: 'category'
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

  public onDgeaVsAllChange(): void {
    this.initDgeaSelection();
    setTimeout(() => this.renderDgeaHeatmap(), 0);
  }

  // Returns all available cell type levels for the dropdown selectors
  getSelectedDgeaComparison(): any | null {
    if (!this.selectedDgeaGroup1) return null;

    const comps = this.getDgeaComparisonMap();

    if (this.dgeaVsAll) {
      const vsAllId = this.makeComparisonId(this.selectedDgeaGroup1, 'all');
      return comps[vsAllId] ?? null;
    }

    if (!this.selectedDgeaGroup2) return null;
    if (this.selectedDgeaGroup1 === this.selectedDgeaGroup2) return null;

    const directId = this.makeComparisonId(this.selectedDgeaGroup1, this.selectedDgeaGroup2);
    if (comps[directId]) {
      return comps[directId];
    }

    const reverseId = this.makeComparisonId(this.selectedDgeaGroup2, this.selectedDgeaGroup1);
    const reverseCmp = comps[reverseId];
    if (!reverseCmp) return null;

    return this.flipDgeaComparison(reverseCmp);
  }

  private flipDgeaComparison(cmp: any): any {
    const flippedTable: any = {};

    if (cmp?.table) {
      for (const col of Object.keys(cmp.table)) {
        if (col === 'avg_log2FC' || col === 'avg_logFC') {
          flippedTable[col] = {};
          for (const gene of Object.keys(cmp.table[col])) {
            const val = Number(cmp.table[col][gene]);
            flippedTable[col][gene] = Number.isFinite(val) ? -val : cmp.table[col][gene];
          }
        } else if (col === 'pct1') {
          flippedTable['pct1'] = cmp.table['pct2'] ? { ...cmp.table['pct2'] } : {};
        } else if (col === 'pct2') {
          flippedTable['pct2'] = cmp.table['pct1'] ? { ...cmp.table['pct1'] } : {};
        } else {
          flippedTable[col] = { ...cmp.table[col] };
        }
      }
    }

    const flippedHeatmap = cmp?.heatmap_context
      ? {
        ...cmp.heatmap_context,
        groups: Array.isArray(cmp.heatmap_context.groups)
          ? [
            this.selectedDgeaGroup1,
            this.selectedDgeaGroup2,
            ...cmp.heatmap_context.groups.filter(
              (g: string) => g !== this.selectedDgeaGroup1 && g !== this.selectedDgeaGroup2
            )
          ]
          : cmp.heatmap_context.groups
      }
      : null;

    return {
      ...cmp,
      group1: this.selectedDgeaGroup1,
      group2: this.selectedDgeaGroup2,
      n1: cmp.n2,
      n2: cmp.n1,
      name: `${this.selectedDgeaGroup1} vs ${this.selectedDgeaGroup2}`,
      table: flippedTable,
      heatmap_context: flippedHeatmap
    };
  }

  // Initialize default selections for the DGEA comparison dropdowns
  initDgeaSelection(): void {
    const levels = this.getDgeaLevels();
    if (!levels.length) return;

    if (!this.selectedDgeaGroup1) {
      this.selectedDgeaGroup1 = levels[0];
    }

    if (this.dgeaVsAll) {
      this.selectedDgeaGroup2 = null;
      return;
    }

    if (!this.selectedDgeaGroup2) {
      const firstDifferent = levels.find(x => x !== this.selectedDgeaGroup1);
      this.selectedDgeaGroup2 = firstDifferent ?? null;
    }
  }


  private loadAndRenderData(dataPath: string, compare: boolean = false): void {

    const tokenType = compare ? 'hexagons_compare' : 'hexagons_main';
    const token = this.nextRequestToken(tokenType);

    // Load GeoJSON data - prepend base URL if it starts with /api/
    const fullUrl = dataPath.startsWith('/api/')
      ? `${this.sessionService.apiUrl}${dataPath}`
      : dataPath;

    console.log('[DEBUG] Loading GeoJSON from:', fullUrl, 'compare=', compare);

    // Use fetch instead of d3.json to handle credentials properly
    fetch(fullUrl, { credentials: 'include' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`HTTP Error: ${response.status} ${response.statusText}`);
        }
        return response.json();
      })
      .then(async (data: GeoJsonData) => {
        if (token !== this.requestTokens[tokenType]) {
          console.warn(`[Race condition prevented] Ignoring stale GeoJSON fetch for ${dataPath}`);
          return;
        }
        console.log('Data loaded from:', fullUrl, data);

        if (!data) {
          throw new Error('Failed to load GeoJSON data');
        }

        compare ? this.geoDataService.setCompareData(data) : this.geoDataService.setData(data);

        const datasetId = compare ? this.selectedDatasetCompare?.id : this.selectedDataset?.id;
        if (compare) {
          this.regulatoryObsmKeysCompare = await this.fetchRegulatoryObsmKeys(datasetId);
        } else {
          this.regulatoryObsmKeysMain = await this.fetchRegulatoryObsmKeys(datasetId);
        }

        console.log('[Regulatory availability] main obsm keys:', this.regulatoryObsmKeysMain);

        // Infer properties for coloring dropdown
        if (compare) {
          this.createHexagonPlot('#hexbin-compare');
          this.compareFeatures = data.features;
          if (this.g_compare && this.currentPathGenerator) {
            this.g_compare
              .style('cursor', 'pointer')
              .style('pointer-events', null)
              .selectAll<SVGPathElement, CellFeature>('path')
              .data(this.compareFeatures, (d: any) => d.properties.barcode)
              .join('path')
              .attr('d', (d: CellFeature) => this.currentPathGenerator!(d) || '')
              .attr('fill', (d: CellFeature) => {
                const value = this.leidenCentralityProps.includes(this.selectedCompareView)
                  ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[this.selectedCompareView]
                  : d.properties?.[this.selectedCompareView];
                if (this.currentCompareLegendType === 'categorical') {
                  return this.colorScaleCompare(String(value));
                } else {
                  const num = this.toNumber(value);
                  return Number.isFinite(num)
                    ? this.continuousColorScaleCompare(num)
                    : '#ccc';
                }
              })
              .style('opacity', 0.8)
              .on('mouseover', (event, d) => this.mouseOver(event, d))
              .on('mouseleave', (event, d) => this.mouseLeave(event, d))
              .on('click', (event, d) => this.displayCellDetails(event, d, true));
          }

        } else {
          this.features = data.features;
        }

        if (data.meta) {
          // save all features for potential downsampling in Xenium datasets
          if (!compare) {
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

            if (compare) {
              this.metaCompare = data.meta;
              this.updatePropertyAvailability('compare');
            } else {
              this.meta = data.meta;
              this.updatePropertyAvailability('main');
            }

            const availableObsCols = this.getDgeaObsCols();
            if (availableObsCols.length && !availableObsCols.includes(this.selectedDgeaObsCol)) {
              this.selectedDgeaObsCol = availableObsCols[0];
            }
            this.dgeaReady = !!this.meta?.['dgea']?.['cell_type'];
            if (this.dgeaReady) {
              this.initDgeaSelection();
            }

            const leidenClusterAnnotations = compare ? this.metaCompare?.['leiden_cluster_annotations'] : this.meta?.['leiden_cluster_annotations'];
            if (leidenClusterAnnotations && typeof leidenClusterAnnotations === 'object') {
              this.clusterCount = Object.keys(leidenClusterAnnotations).length;
            }

            const interval = compare ? this.metaCompare?.['interval'] : this.meta?.['interval'];
            if (Array.isArray(interval) && interval.length > 0) {
              this.maxInterval = interval.length - 1;
            }
            this.selectedRegulatoryScore =
              compare ? this.metaCompare['grn_score_names']?.[0] : this.meta['grn_score_names']?.[0] || null;
            if (compare) {
              this.currentCompareLegendType = this.isContinuousScale(this.selectedCompareView, this.compareFeatures, true) ? 'continuous' : 'categorical';
              this.geneSetsGenie3Compare = this.metaCompare['genie_genesets'] || {};
              this.geneSetsSpongeCompare = this.metaCompare['sponge_genesets'] || {};
            } else {
              this.currentLegendType = this.isContinuousScale(this.selectedView, this.features, false) ? 'continuous' : 'categorical';
              this.geneSetsGenie3 = this.meta['genie_genesets'] || {};
              this.geneSetsSponge = this.meta['sponge_genesets'] || {};
            }
            // Populate dropdown options from gene set keys
            this.genie3Elements = Object.keys(this.geneSetsGenie3);
            this.spongeElements = Object.keys(this.geneSetsSponge);
            this.selectedGeneSetGenie3 =
              Object.keys(compare ? this.metaCompare['genie_genesets'] || {} : this.meta['genie_genesets'] || {})[0] || null;
            this.selectedGeneSetSponge =
              Object.keys(compare ? this.metaCompare['sponge_genesets'] || {} : this.meta['sponge_genesets'] || {})[0] || null;

            if (compare) {
              this.isLoadingGenie3Compare = !!this.selectedGeneSetGenie3Compare;
              this.isLoadingSpongeCompare = !!this.selectedGeneSetSpongeCompare;
              this.previousGeneSetGenie3Compare = null;
              this.previousGeneSetSpongeCompare = null;
              this.selectedRegulatoryScoreCompare = this.metaCompare['grn_score_names']?.[0] || null;
            } else {
              this.isLoadingGenie3 = !!this.selectedGeneSetGenie3;
              this.isLoadingSponge = !!this.selectedGeneSetSponge;
              this.previousGeneSetGenie3 = null;
              this.previousGeneSetSponge = null;
              this.selectedRegulatoryScore = this.meta['grn_score_names']?.[0] || null;
            }
            // Trigger initial graph visualization
            setTimeout(() => {
              this.onGeneSetChange(compare);
            }, 100);
          }
        }
        const featuresForProps = compare ? this.compareFeatures : this.features;
        const firstProps = featuresForProps?.[0]?.properties || {};


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
        compare ? this.colorablePropertiesCompare = Array.from(allKeys) : this.colorableProperties = Array.from(allKeys);
        console.log('Initial colourable properties (including potential ones):',
          compare ? this.colorablePropertiesCompare : this.colorableProperties);

        // Alphabetical order
        compare ? this.colorablePropertiesCompare.sort((a, b) => a.localeCompare(b)) : this.colorableProperties.sort((a, b) => a.localeCompare(b));

        // Group similar properties together
        const chromvarKeys = ['chromvar_total_sum'];
        const currentgroupedProperties = [
          { key: 'Scores', value: this.colorableProperties.filter((p) => scoreKeys.includes(p)) },
          { key: 'LIANA+', value: this.colorableProperties.filter((p) => lianaKeys.includes(p)) },
          { key: 'ChromVAR', value: this.colorableProperties.filter((p) => chromvarKeys.includes(p)) },
          {
            key: 'Other', value: this.colorableProperties.filter(
              (p) => !scoreKeys.includes(p) && !lianaKeys.includes(p) && !chromvarKeys.includes(p)
            )
          },
        ];

        compare ? this.groupedPropertiesCompare = currentgroupedProperties : this.groupedProperties = currentgroupedProperties;


        // compute availability for everything we've decided to show
        const candidates = new Set(compare ? this.colorablePropertiesCompare : this.colorableProperties);

        const targetAvailabilityMap: { [prop: string]: boolean } = {};

        candidates.forEach((prop) => {
          if (prop === 'regulatory_scores') {
            targetAvailabilityMap[prop] = this.hasRegulatoryScoresData(compare ? 'compare' : 'main');
          } else if (prop === 'gene_expression') {
            targetAvailabilityMap[prop] = this.hasGeneExpressionData(compare ? 'compare' : 'main');
          } else if (this.leidenCentralityProps.includes(prop)) {
            targetAvailabilityMap[prop] = (compare ? this.compareFeatures : this.features).some((f) => {
              // FIX: Added 'compare' argument to getLeidenClusterAnnotation
              const val = this.getLeidenClusterAnnotation(f.properties.leiden, compare)?.centrality?.[prop];
              return val !== undefined && val !== null && val !== '';
            });
          } else {
            targetAvailabilityMap[prop] = (compare ? this.compareFeatures : this.features).some((f) => {
              const val = f.properties ? f.properties[prop] : undefined;
              return val !== undefined && val !== null && val !== '';
            });
          }
        });

        // Conditionally assign the calculated map to the correct view
        if (compare) {
          this.propertyAvailabilityCompare = targetAvailabilityMap;
        } else {
          this.propertyAvailability = targetAvailabilityMap;
        }
        console.log(`Property availability map (${compare ? 'compare' : 'main'}):`, targetAvailabilityMap);

        // rebuild grouped properties with full property list; template handles disabled state
        this.groupedProperties = [
          { key: 'Scores', value: this.colorableProperties.filter((p) => scoreKeys.includes(p)) },
          { key: 'LIANA+', value: this.colorableProperties.filter((p) => lianaKeys.includes(p)) },
          { key: 'ChromVAR', value: this.colorableProperties.filter((p) => chromvarKeys.includes(p)) },
          {
            key: 'Other', value: this.colorableProperties.filter(
              (p) => !scoreKeys.includes(p) && !lianaKeys.includes(p) && !chromvarKeys.includes(p)
            ),
          },
        ];

        const firstAvailableProperty = compare
          ? this.colorablePropertiesCompare.find((p) => this.propertyAvailable(p, 'compare'))
          : this.colorableProperties.find((p) => this.propertyAvailable(p, 'main'));

        if (compare) {
          if (this.propertyAvailable('regulatory_scores', 'compare')) {
            this.selectedCompareView = 'regulatory_scores';
          } else if (this.propertyAvailable('cell_type', 'compare')) {
            this.selectedCompareView = 'cell_type';
          } else if (firstAvailableProperty) {
            this.selectedCompareView = firstAvailableProperty;
          } else {
            // Keep empty when nothing is available.
            this.selectedCompareView = this.colorablePropertiesCompare[0] || '';
          }
          this.currentCompareLegendType = this.isContinuousScale(this.selectedCompareView, this.compareFeatures, true) ? 'continuous' : 'categorical';
        } else {

          if (this.propertyAvailable('regulatory_scores', 'main')) {
            this.selectedView = 'regulatory_scores';
          } else if (this.propertyAvailable('cell_type', 'main')) {
            this.selectedView = 'cell_type';
          } else if (firstAvailableProperty) {
            this.selectedView = firstAvailableProperty;
          } else {
            // Keep empty when nothing is available.
            this.selectedView = this.colorableProperties[0] || '';
          }
          this.currentLegendType = this.isContinuousScale(this.selectedView, this.features, false) ? 'continuous' : 'categorical';
        }



        const width = 1200;
        const height = 1000;

        const projection = d3.geoIdentity().fitSize([width, height], {
          type: 'FeatureCollection',
          features: compare ? this.compareFeatures : (this.isXenium ? this.fullFeatures : this.features),
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
              const value = this.leidenCentralityProps.includes(compare ? this.selectedCompareView : this.selectedView)
                ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[compare ? this.selectedCompareView : this.selectedView]
                : d.properties?.[compare ? this.selectedCompareView : this.selectedView];
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
          (compare ? this.g_compare : this.g)
            .style('cursor', 'pointer')
            .style('pointer-events', null)
            .selectAll<SVGPathElement, CellFeature>('path')
            .data(compare ? this.compareFeatures : this.features, (d: any) => d.properties.barcode)
            .join('path')
            .attr('d', (d: CellFeature) => pathGenerator(d) || '')
            .attr('fill', (d: CellFeature) => {
              const value = this.leidenCentralityProps.includes(compare ? this.selectedCompareView : this.selectedView)
                ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[compare ? this.selectedCompareView : this.selectedView]
                : d.properties?.[compare ? this.selectedCompareView : this.selectedView];
              if (compare ? this.currentCompareLegendType === 'categorical' : this.currentLegendType === 'categorical') {
                return this.colorScale(String(value));
              } else {
                const num = this.toNumber(value);
                return Number.isFinite(num)
                  ? compare ? this.continuousColorScaleCompare(num)
                    : this.continuousColorScale(num)
                  : '#ccc';
              }
            })
            .style('opacity', 0.8)
            .on('mouseover', (event, d) => this.mouseOver(event, d))
            .on('mouseleave', (event, d) => this.mouseLeave(event, d))
            .on('click', (event, d) => this.displayCellDetails(event, d, compare));
        }

        this.onColorbyPropertyChange(compare);
        setTimeout(() => {
          if (compare) {
            this.isLoadingCompare = false;
          } else {
            this.isLoadingHexagons = false;
          }
          this.checkInitializationComplete(compare);
        }, 0);
      })
      .catch((error) => {
        console.error('Error loading or rendering data:', error);
        if (compare) {
          this.isLoadingCompare = false;
        } else {
          this.isLoadingHexagons = false;
        }
        this.checkInitializationComplete(compare);
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

  private normalizeGeneKey(gene: string): string {
    return gene.trim().toLowerCase();
  }

  private buildGeneDomainCacheKey(datasetId: string, gene: string): string {
    return `${datasetId}::${this.normalizeGeneKey(gene)}`;
  }

  private getGeneDomainRequests(): Array<{ datasetId: string; gene: string; cacheKey: string }> {
    const requests: Array<{ datasetId: string; gene: string; cacheKey: string }> = [];
    const seen = new Set<string>();

    const add = (datasetId: string | undefined, gene: string | null) => {
      const cleanGene = (gene ?? '').trim();
      if (!datasetId || !cleanGene) return;
      const cacheKey = this.buildGeneDomainCacheKey(datasetId, cleanGene);
      if (seen.has(cacheKey)) return;
      seen.add(cacheKey);
      requests.push({ datasetId, gene: cleanGene, cacheKey });
    };

    if (this.selectedView === 'gene_expression') {
      add(this.selectedDataset?.id, this.selectedGeneExpressionMain);
    }

    if (this.compareMode && this.selectedCompareView === 'gene_expression') {
      const compareDatasetId = this.selectedDatasetCompare?.id ?? this.selectedDataset?.id;
      const compareGene = this.selectedGeneExpressionCompare ?? this.selectedGeneExpressionMain;
      add(compareDatasetId, compareGene);
    }

    return requests;
  }

  private getGeneDomainContextKey(requests: Array<{ cacheKey: string }>): string {
    return requests.map(r => r.cacheKey).sort().join('|');
  }

  private getSharedDomainForGeneExpressionView(): { min: number; max: number } | null {
    return this.sharedGeneExpressionDomain;
  }

  private cloneGeoJsonData(data: GeoJsonData | null): GeoJsonData | null {
    if (!data) return null;

    try {
      if (typeof structuredClone === 'function') {
        return structuredClone(data);
      }
    } catch (e) {
      console.warn('[compare:data] structuredClone failed, falling back to JSON clone', e);
    }

    return JSON.parse(JSON.stringify(data)) as GeoJsonData;
  }

  private async loadCompareGeoJsonData(): Promise<GeoJsonData | null> {
    const comparePath = this.selectedDatasetCompare?.geojson_path
      ?? this.selectedDataset?.geojson_path
      ?? this.dataPath;

    if (!comparePath) {
      console.warn('[compare init] no compare geojson path available');
      return null;
    }

    const fullUrl = comparePath.startsWith('/api/')
      ? `${this.sessionService.apiUrl}${comparePath}`
      : comparePath;

    const response = await fetch(fullUrl, { credentials: 'include' });
    if (!response.ok) {
      throw new Error(`HTTP Error: ${response.status} ${response.statusText} (${fullUrl})`);
    }

    const data = await response.json() as GeoJsonData;
    console.log('[compare init] loaded geojson from path', {
      comparePath,
      fullUrl,
      featureCount: data?.features?.length ?? 0,
      selectedCompareDatasetId: this.selectedDatasetCompare?.id ?? null,
      selectedMainDatasetId: this.selectedDataset?.id ?? null,
    });

    return data;
  }

  private extractViewValue(feature: CellFeature, view: string): unknown {
    if (this.leidenCentralityProps.includes(view)) {
      return this.getLeidenClusterAnnotation(feature.properties.leiden)?.centrality?.[view];
    }
    return feature.properties[view];
  }

  private collectFiniteValuesForView(features: CellFeature[], view: string): number[] {
    return features
      .map((f) => this.toNumber(this.extractViewValue(f, view)))
      .filter((n) => Number.isFinite(n));
  }

  private getMinMaxForView(features: CellFeature[], view: string): { min: number; max: number } | null {
    const values = this.collectFiniteValuesForView(features, view);
    if (!values.length) {
      console.log('[domain:getMinMaxForView] no finite values', {
        view,
        featureCount: features.length,
      });
      return null;
    }

    const result = {
      min: Math.min(...values),
      max: Math.max(...values),
    };

    console.log('[domain:getMinMaxForView] computed', {
      view,
      featureCount: features.length,
      finiteCount: values.length,
      min: result.min,
      max: result.max,
    });

    return result;
  }

  private getPairedContinuousDomainForCompare(): { min: number; max: number } | null {
    if (!this.compareMode) {
      console.log('[domain:getPaired] skipped: compare mode is off');
      return null;
    }

    const mainView = this.selectedView;
    const compareView = this.selectedCompareView;
    const mainFeatures = this.features || [];
    const compareFeatures = this.compareFeatures || [];
    if (!mainFeatures.length || !compareFeatures.length) return null;

    const mainIsContinuous = this.isContinuousScale(mainView, mainFeatures);
    const compareIsContinuous = this.isContinuousScale(compareView, compareFeatures);
    if (!mainIsContinuous || !compareIsContinuous) {
      console.log('[domain:getPaired] skipped: one side is not continuous', {
        mainView,
        compareView,
        mainIsContinuous,
        compareIsContinuous,
      });
      return null;
    }

    const mainMinMax = this.getMinMaxForView(mainFeatures, mainView);
    const compareMinMax = this.getMinMaxForView(compareFeatures, compareView);
    if (!mainMinMax || !compareMinMax) {
      console.log('[domain:getPaired] skipped: missing min/max', {
        mainView,
        compareView,
        mainMinMax,
        compareMinMax,
      });
      return null;
    }

    const merged = {
      min: Math.min(mainMinMax.min, compareMinMax.min),
      max: Math.max(mainMinMax.max, compareMinMax.max),
    };

    console.log('[domain:getPaired] merged domain', {
      mainView,
      compareView,
      mainMinMax,
      compareMinMax,
      merged,
    });

    return merged;
  }

  private shouldUseComparablePairedDomain(): boolean {
    return this.compareMode && this.colorByProperty === this.selectedCompareView;
  }

  private applySharedDomainAndRepaint(domain: { min: number; max: number } | null, contextKey: string, token: number): void {
    if (token !== this.geneDomainToken) return;
    this.sharedGeneExpressionDomain = domain;
    this.sharedGeneExpressionContextKey = domain ? contextKey : null;

    this.updateHexColors();
    if (this.compareMode) {
      this.updateHexColors('#hexbin-compare');
    }
  }

  private async refreshSharedGeneExpressionDomain(): Promise<void> {
    const requests = this.getGeneDomainRequests();
    if (!requests.length) {
      this.sharedGeneExpressionDomain = null;
      this.sharedGeneExpressionContextKey = null;
      return;
    }

    const contextKey = this.getGeneDomainContextKey(requests);
    const token = ++this.geneDomainToken;
    const now = Date.now();

    const pending = requests.filter((r) => {
      const cached = this.geneDomainCache.get(r.cacheKey);
      return !cached || cached.expiresAt <= now;
    });

    if (pending.length > 0) {
      const results = await Promise.allSettled(
        pending.map((r) => firstValueFrom(
          this.sessionService.callWithSession(() =>
            this.http.get<GeneStatsResponse>(
              `${this.sessionService.apiUrl}/X_stats/${encodeURIComponent(r.gene)}?dataset_ids=${encodeURIComponent(r.datasetId)}`,
              { withCredentials: true },
            )
          )
        ))
      );

      if (token !== this.geneDomainToken) {
        return;
      }

      results.forEach((result, idx) => {
        const req = pending[idx];
        if (result.status === 'fulfilled') {
          this.geneDomainCache.set(req.cacheKey, {
            min: result.value.global_min,
            max: result.value.global_max,
            expiresAt: now + this.geneDomainCacheTtlMs,
          });
        } else {
          console.warn('[gene-domain] failed to fetch global stats', req, result.reason);
        }
      });
    }

    const valid = requests
      .map((r) => this.geneDomainCache.get(r.cacheKey))
      .filter((v): v is { min: number; max: number; expiresAt: number } => !!v && v.expiresAt > Date.now());

    if (!valid.length) {
      this.applySharedDomainAndRepaint(null, contextKey, token);
      return;
    }

    const min = Math.min(...valid.map((v) => v.min));
    const max = Math.max(...valid.map((v) => v.max));
    this.applySharedDomainAndRepaint({ min, max }, contextKey, token);
  }

  public onGeneExpressionSelected(event: { gene: string; action: string }, panel: 'main' | 'compare'): void {
    if (event.action !== 'gene_expression') return;
    const selectedGene = String(event.gene ?? '').trim();
    if (!selectedGene) return;

    if (panel === 'main') {
      this.selectedGeneExpressionMain = selectedGene;
    } else {
      this.selectedGeneExpressionCompare = selectedGene;
    }

    this.refreshSharedGeneExpressionDomain();
  }

  public onColorbyPropertyChange(compare: boolean = false): void {
    const colorProp = compare ? this.selectedCompareView : this.selectedView;
    const regulatoryScore = compare ? this.selectedRegulatoryScoreCompare : this.selectedRegulatoryScore;
    const geneSetGenie3 = compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3;
    const geneSetSponge = compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge;

    console.log('[onColorbyPropertyChange] colorByProperty changed to:', colorProp);

    const dualGeneExpressionCompare =
      this.compareMode &&
      this.selectedView === 'gene_expression' &&
      this.selectedCompareView === 'gene_expression';

    if (dualGeneExpressionCompare) {
      this.refreshSharedGeneExpressionDomain();
    } else {
      this.sharedGeneExpressionDomain = null;
      this.sharedGeneExpressionContextKey = null;
    }

    if (colorProp === 'regulatory_scores') {
      if (regulatoryScore?.endsWith('genie3') && geneSetGenie3) {
        this.fetchAndUpdate(regulatoryScore, geneSetGenie3, compare);
        this.updateSubgraphGenie3(compare);
        this.updateSubgraphSponge(compare);
      } else if (regulatoryScore?.endsWith('sponge') && geneSetSponge) {
        this.fetchAndUpdate(regulatoryScore, geneSetSponge, compare);
        this.updateSubgraphGenie3(compare);
        this.updateSubgraphSponge(compare);
      }
    }

    if (compare) {
      console.log('[onColorbyPropertyChange] Updating compare hexagons');
      this.updateHexColors('#hexbin-compare');
    } else {
      console.log('[onColorbyPropertyChange] Updating main hexagons');
      this.updateHexColors();
    }
  }

  /**
   * Determine whether the given view/property should be treated as a continuous numeric scale.
   * If `features` is provided, that dataset is used for the test (used for compare view). Otherwise
   * the main `this.features` is used.
   */
  isContinuousScale(view?: string, features?: CellFeature[], compare: boolean = false): boolean {
    const property = view ?? (compare ? this.selectedCompareView : this.selectedView);
    const sourceFeatures = Array.isArray(features)
      ? features
      : (compare ? this.compareFeatures : this.features) || [];

    const valuesRaw = sourceFeatures.map((f) => {
      if (this.leidenCentralityProps.includes(property)) {
        const clusterAnnotation = this.getLeidenClusterAnnotation(f.properties.leiden);
        return clusterAnnotation?.centrality?.[property];
      }
      return f.properties[property];
    });

    // FILTER OUT undefined/null values so they don't evaluate to NaN and break the check
    const validRaw = valuesRaw.filter(v => v !== undefined && v !== null && v !== '');

    // If completely empty (e.g., initial state), default to continuous
    if (validRaw.length === 0) return true;

    const numericValues = validRaw.map((v) => this.toNumber(v));
    const allNumbers = numericValues.every((n) => Number.isFinite(n));
    const allIntegers = allNumbers && numericValues.every((n) => Number.isInteger(n));
    const uniqueIntegerCount = allIntegers ? new Set(numericValues).size : 0;
    const shouldTreatAsCategorical = allIntegers && uniqueIntegerCount <= 20;

    return !shouldTreatAsCategorical;
  }

  private getViewVariablesToUpdate(containerName: string) {
    const isMainView = containerName === '#hexbin';
    const viewToUpdate = isMainView ? this.selectedView : this.selectedCompareView;
    const featuresToUpdate = isMainView ? this.features : (this.compareFeatures || []);
    const gToUpdate = isMainView ? this.g : this.g_compare;
    const gPathsToUpdate = isMainView ? this.g_paths : this.g_paths_compare;
    const ordinalScaleToUpdate = isMainView ? this.colorScale : this.colorScaleCompare;
    const continuousScaleToUpdate = isMainView ? this.continuousColorScale : this.continuousColorScaleCompare;
    // Pass the compare flag to isContinuousScale for correct logic
    const isContinuous = this.isContinuousScale(viewToUpdate, featuresToUpdate, !isMainView);
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

  // This function checks if the selected view/property is available in the features. If not, it falls back to a default property or the first available property.

  private getAvailableView(selectedView: string, features: CellFeature[], fallback: string = 'cell_type'): string {

    if (!features.length) return fallback;

    // 1. Always trust dynamically fetched properties (even if sparse/empty on cell 0)
    const dynamicProps = [
      'gene_expression', 'regulatory_scores', 'ligand_receptor_relationships',
      'cell_comp_tf_activity_similarity', 'tf_activity', 'pathway_activity', 'chromvar_total_sum'
    ];
    if (dynamicProps.includes(selectedView)) {
      return selectedView;
    }

    // 2. For static properties, check if it exists on the first feature
    const firstFeature = features[0];
    const hasProperty = this.leidenCentralityProps.includes(selectedView)
      ? true
      : selectedView in (firstFeature.properties || {});

    if (hasProperty) return selectedView;
    if (fallback in (firstFeature.properties || {})) return fallback;

    const keys = Object.keys(firstFeature.properties || {});
    return keys.length ? keys[0] : fallback;
  }

  public updateHexColors(containerName?: string): void {
    if (!containerName) {
      containerName = '#hexbin';
    }

    // Determine which features and selected view to use
    const isMainView = containerName === '#hexbin';
    const features = isMainView ? this.features : this.compareFeatures;
    const selectedView = isMainView ? this.selectedView : this.selectedCompareView;

    // Use getAvailableView to ensure the property exists
    const viewToUse = this.getAvailableView(selectedView, features, 'cell_type');
    if (viewToUse !== selectedView) {
      console.warn(`[updateHexColors] Property "${selectedView}" not found, falling back to "${viewToUse}".`);
      // Optionally, update the selectedView variable here if you want the UI to reflect the fallback
      // if (isMainView) this.selectedView = viewToUse; else this.selectedCompareView = viewToUse;
    }

    // Get the rest of the view variables (scales, g, etc.)
    const viewVariablesToUpdate = this.getViewVariablesToUpdate(containerName);

    // Overwrite the view property for this rendering
    // (since viewVariablesToUpdate.view is readonly, just use viewToUse in all logic below)
    let valuesRaw: any[] = [];
    if (features.length > 0) {
      const firstFeature = features[0];
      const hasProperty = this.leidenCentralityProps.includes(viewToUse)
        ? true
        : viewToUse in (firstFeature.properties || {});
      console.log('[updateHexColors] First feature properties keys:', Object.keys(firstFeature.properties || {}));
      console.log('[updateHexColors] Looking for property:', viewToUse, '- Exists:', hasProperty);
      if (!hasProperty) {
        console.warn('[updateHexColors] Property not found in features! Available:', Object.keys(firstFeature.properties || {}));
        // Optionally, show a user-facing message here
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
      .data(features);

    if (containerName === '#hexbin-compare' && this.compareMode && sel.size() === 0) {
      console.warn('[updateHexColors] Compare layer has no paths; re-initializing compare plot');
      this.loadAndRenderData(this.selectedDatasetCompare?.geojson_path || '', true);
      return;
    }

    valuesRaw = features.map((f) => this.extractViewValue(f, viewToUse));
    const numericValues = valuesRaw.map((v) => this.toNumber(v));

    // Determine if the property is continuous
    const isContinuous = this.isContinuousScale(viewToUse, features, !isMainView);

    if (isContinuous) {
      const compareSharedDomain = this.getPairedContinuousDomainForCompare();
      console.log('[updateHexColors] Computed compareSharedDomain:', compareSharedDomain);
      const sharedDomain = compareSharedDomain ?? (viewToUse === 'gene_expression'
        ? this.getSharedDomainForGeneExpressionView()
        : null;

      const domainSource = sharedDomain ? 'shared-gene-expression-domain' : 'local-current-view-domain';

      let min = sharedDomain ? sharedDomain.min : Math.min(...numericValues);
      let max = sharedDomain ? sharedDomain.max : Math.max(...numericValues);

      if (!Number.isFinite(min) || !Number.isFinite(max)) {
        min = 0;
        max = 1;
      }

      console.log('[updateHexColors] applying continuous domain', {
        container: containerName,
        view: viewVariablesToUpdate.view,
        compareMode: this.compareMode,
        mainView: this.colorByProperty,
        compareView: this.selectedCompareView,
        domainSource,
        localFiniteValueCount: numericValues.length,
        min,
        max,
      });
      if (min === max) {
        const eps = min === 0 ? 1 : Math.abs(min) * 0.01;
        min -= eps;
        max += eps;
      }

      // Always normalize to [0, 1] for consistent color usage across all views
      this.continuousColorScale.domain([0, 1]);
      this.continuousColorScaleCompare.domain([0, 1]);
      viewVariablesToUpdate.continuous.domain([0, 1]);

      // Legend always shows [0, 1]
      this.currentLegendDomain = [0, 1];
      this.currentLegendDomainCompare = [0, 1];
      this.currentLegendType = 'continuous';
      this.currentCompareLegendType = 'continuous';
      viewVariablesToUpdate.setLegendDomain([0, 1]);
      viewVariablesToUpdate.setLegendType('continuous');

      sel
        .transition()
        .duration(300)
        .attr('stroke-width', 1)
        .attr('stroke', 'transparent')
        .attr('fill', (d) => {
          const raw = this.extractViewValue(d, viewToUse);
          const n = this.toNumber(raw);
          return Number.isFinite(n)
            ? viewVariablesToUpdate.continuous(n)
            : '#ccc';
        });
    } else {
      // categorical scale
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
          const raw = this.extractViewValue(d, viewToUse);
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

  public updateSubgraphGenie3(compare: boolean = false): void {
    const token = this.nextRequestToken('genie3');
    console.log('Updating AUCELL graph for Genie3...');
    compare ? this.isLoadingGenie3 = true : this.isLoadingGenie3Compare = true;
    const geneSet = compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3;
    const datasetId = compare ? this.selectedDatasetCompare?.id : this.selectedDataset?.id;
    const graphContainerId = compare ? '#aucell_graph_genie3_compare' : '#aucell_graph_genie3';
    d3.select(graphContainerId).selectAll('*').remove();

    if (!geneSet || !datasetId) {
      compare ? this.isLoadingGenie3 = false : this.isLoadingGenie3Compare = false;
      return;
    }

    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/geneset_connections_genie?gene_set_name=${encodeURIComponent(geneSet || '')}&dataset_id=${encodeURIComponent(datasetId || '')}`,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          if (token !== this.requestTokens['genie3']) return;
          {
            const payload: any = res;
            const data = payload['connections'] as { regulatoryGene: string; targetGene: string; weight: number }[];
            compare ? this.genie3NetworkCompare = data.map((d) => ({
              source: d.regulatoryGene,
              target: d.targetGene,
              weight: d.weight,
            })) : this.genie3Network = data.map((d) => ({
              source: d.regulatoryGene,
              target: d.targetGene,
              weight: d.weight,
            }));

            // slider_data may be returned as an object or an array; handle both safely
            const sliderData: any = payload['slider_data'];

            if (sliderData && typeof sliderData === 'object') {
              if (compare) {
                this.genie3SliderDataCompare = {
                  step: sliderData.step || 1,
                  min_border: sliderData.min_border || 0,
                  max_border: sliderData.max_border || 100,
                  default_value: sliderData.default_value || 50,
                };
                this.genie3WeightCutoffCompare = this.genie3SliderDataCompare.default_value;
              }
              else {
                this.genie3SliderData = {
                  step: sliderData.step || 1,
                  min_border: sliderData.min_border || 0,
                  max_border: sliderData.max_border || 100,
                  default_value: sliderData.default_value || 50,
                };
                this.genie3WeightCutoff = this.genie3SliderData.default_value;
              }
            }

            this.visualizeGenie3Subgraph(compare);

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

  public visualizeGenie3Subgraph(compare: boolean = false): void {
    const graphContainerId = compare ? '#aucell_graph_genie3_compare' : '#aucell_graph_genie3';
    const geneSet = compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3;


    compare ? this.isLoadingGenie3Compare = true : this.isLoadingGenie3 = true;
    d3.select(graphContainerId).html('');

    if (!geneSet || compare ? !this.genie3NetworkCompare : !this.genie3Network || compare ? this.genie3NetworkCompare.length === 0 : this.genie3Network.length === 0) {
      compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
      this.checkInitializationComplete(compare);
      return;
    }

    let regulator = geneSet ?? '';
    let targets = this.geneSetsGenie3[regulator] || [];

    let nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    let edges: { source: string; target: string; weight: number }[] = [];

    let candidateEdges: { source: string; target: string; weight: number }[] = [];
    let slicedEdges: { source: string; target: string; weight: number }[] = [];

    if (compare) {
      candidateEdges = this.genie3NetworkCompare.filter((edge) => edge.weight > this.genie3WeightCutoffCompare).map((e) => ({
        source: String(e.source),
        target: String(e.target),
        weight: e.weight,
      }));
    }
    else {
      candidateEdges = this.genie3Network.filter((edge) => edge.weight > this.genie3WeightCutoff).map((e) => ({
        source: String(e.source),
        target: String(e.target),
        weight: e.weight,
      }));
    }
    candidateEdges.sort((a, b) => b.weight - a.weight);
    slicedEdges = candidateEdges.slice(0, this.genie3MinEdges);

    // Infer nodes from edges
    const nodeSet = new Set<string>();

    slicedEdges.forEach((edge) => {
      if (!(nodeSet.has(edge.source) && nodeSet.has(edge.target))) {
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
      .select(graphContainerId)
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

    setTimeout(() => this.updateGraphWidths(), 50);

    compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
    this.checkInitializationComplete(compare);


  }


  public updateSubgraphSponge(compare: boolean = false): void {
    const graphContainerId = compare ? '#aucell_graph_sponge_compare' : '#aucell_graph_sponge';
    const datasetId = compare ? this.selectedDatasetCompare?.id : this.selectedDataset?.id;
    const geneSet = compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge;

    const token = this.nextRequestToken('sponge');
    console.log('Updating AUCELL graph for Sponge...');
    d3.select(graphContainerId).selectAll('*').remove();

    if (!geneSet || !datasetId) {
      return;
    }

    compare ? this.isLoadingSpongeCompare = true : this.isLoadingSponge = true;
    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/geneset_connections_sponge?gene_set_name=${encodeURIComponent(geneSet || '')}&dataset_id=${encodeURIComponent(datasetId || '')}`,
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

            this.visualizeSpongeSubgraph(compare);
          }
        },
        error: (err) => {
          if (this.requestTokens['sponge'] !== this.requestTokens['sponge']) {
            console.error(
              `[Backend] Failed to load Sponge Connections for["${geneSet}"]`,
              err,
            );
          }
        },
      });
  }

  public visualizeSpongeSubgraph(compare: boolean = false): void {
    const datasetId = compare ? this.selectedDatasetCompare?.id : this.selectedDataset?.id;
    const geneSet = compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge;
    const graphContainerId = compare ? '#aucell_graph_sponge_compare' : '#aucell_graph_sponge';

    compare ? this.isLoadingSpongeCompare = true : this.isLoadingSponge = true;
    d3.select(graphContainerId).html('');

    if (!geneSet || !datasetId || compare ? !this.spongeNetworkCompare || this.spongeNetworkCompare.length === 0 : !this.spongeNetwork || this.spongeNetwork.length === 0) {
      compare ? this.isLoadingSpongeCompare = false : this.isLoadingSponge = false;
      this.checkInitializationComplete(compare);
      return;
    }

    let regulator = geneSet ?? '';
    let targets = this.geneSetsGenie3[regulator] || [];


    let nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    let edges: { source: string; target: string; p_adjusted: number }[] = [];

    let candidateEdges: { source: string; target: string; p_adjusted: number }[] = [];
    let slicedEdges: { source: string; target: string; p_adjusted: number }[] = [];

    if (compare) {
      candidateEdges = this.spongeNetworkCompare.filter((edge) => edge.p_adjusted < this.spongePValueCutoff).map((e) => ({
        source: String(e.source),
        target: String(e.target),
        p_adjusted: e.p_adjusted,
      }));
    }
    else {
      candidateEdges = this.spongeNetwork.filter((edge) => edge.p_adjusted < this.spongePValueCutoff).map((e) => ({
        source: String(e.source),
        target: String(e.target),
        p_adjusted: e.p_adjusted,
      }));
    }

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
      .select(graphContainerId)
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

    compare ? this.isLoadingSpongeCompare = false : this.isLoadingSponge = false;
    this.checkInitializationComplete(compare);

    setTimeout(() => this.updateGraphWidths(), 50);

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
        (this.selectedView === 'leiden' &&
          d.properties.leiden === this.selectedCell.properties.leiden))
    )
      return;
    if (
      this.selectedCellCompare &&
      d.properties.barcode === this.selectedCellCompare.properties.barcode
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

  public displayCellDetails(event: MouseEvent, cell: CellFeature, compare: boolean = false): void {
    if (!compare) {
      this.resetClusterExtension();
      this.selectedCell = cell;
      if (this.selectedView === 'regulatory_scores') {
        this.getRegulatoryScoresforSpots(cell.properties.barcode, this.selectedDataset?.id)
      }
      if (this.selectedView === 'leiden') {
        this.displayClusterDetails(cell.properties.leiden);
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
    else {
      this.selectedCellCompare = cell;
      if (this.selectedCompareView === 'regulatory_scores') {
        this.getRegulatoryScoresforSpots(
          cell.properties.barcode,
          this.selectedDatasetCompare?.id ?? this.selectedDataset?.id,
        )
      }
      d3.select(event.target as SVGElement)
        .transition()
        .attr('stroke', 'black');
    }
  }


  public displayClusterDetails(clusterId: number): void {
    this.selectedCluster = clusterId;
    this.clusterCells = this.features.filter(
      (cell) => cell.properties.leiden === clusterId,
    );
    this.calculateClusterStats();

    // Initialize co-occurrence table for this cluster
    this.updateCoOccurrenceTable(isCompare);

    const targetClusterCells = isCompare ? this.clusterCellsCompare : this.clusterCells;
    if (targetClusterCells.length > 0) {
      if (isCompare) {
        this.selectedCellCompare = targetClusterCells[0];
      } else {
        this.selectedCell = targetClusterCells[0];
        // Neighborhood enrichment will render when the Cluster Information tab is viewed
        // via onTabChange handler
        setTimeout(() => this.updateSubgraphGenie3(), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
      }
    }
  }

  public selectCluster(clusterId: number, compare: boolean = false): void {
    if (compare) {
      this.selectedCluster = clusterId;
      this.clusterCells = this.compareFeatures.filter(
        (cell) => cell.properties.leiden === clusterId,
      );
      this.calculateClusterStats(); // If this uses clusterCells, it's fine
      this.updateCoOccurrenceTable(); // If this uses clusterCells, it's fine

      if (this.clusterCells.length > 0) {
        this.selectedCellCompare = this.clusterCells[0];
        setTimeout(() => this.updateSubgraphGenie3(true), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDatasetCompare), 100);
      }
    } else {
      this.selectedCluster = clusterId;
      this.clusterCells = this.features.filter(
        (cell) => cell.properties.leiden === clusterId,
      );
      this.calculateClusterStats();
      this.updateCoOccurrenceTable();

      if (this.clusterCells.length > 0) {
        this.selectedCell = this.clusterCells[0];
        setTimeout(() => this.updateSubgraphGenie3(false), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
      }
    }
  }

  public onGeneSetChange(compare: boolean = false): void {
    const geneSetGenie3Changed = compare ? this.selectedGeneSetGenie3 : this.selectedGeneSetGenie3Compare;
    const geneSetSpongeChanged = compare ? this.selectedGeneSetSponge : this.selectedGeneSetSpongeCompare;
    const previousGeneSetGenie3Changed = compare ? this.previousGeneSetGenie3 : this.previousGeneSetGenie3Compare;
    const previousGeneSetSpongeChanged = compare ? this.previousGeneSetSponge : this.previousGeneSetSpongeCompare;
    const regulatoryScoreChanged = compare ? this.selectedRegulatoryScoreCompare : this.selectedRegulatoryScore;
    const graphIdGenie3 = compare ? '#aucell_graph_genie3_compare' : '#aucell_graph_genie3';
    const graphIdSponge = compare ? '#aucell_graph_sponge_compare' : '#aucell_graph_sponge';

    if (geneSetGenie3Changed !== previousGeneSetGenie3Changed) {
      compare ? this.previousGeneSetGenie3Compare = this.selectedGeneSetGenie3 : this.previousGeneSetGenie3 = this.selectedGeneSetGenie3;
      if (geneSetGenie3Changed) {

        d3.select(graphIdGenie3).selectAll('*').remove();
        compare ? this.isLoadingGenie3Compare = true : this.isLoadingGenie3 = true;
        setTimeout(() => {
          this.updateSubgraphGenie3(compare);
          if (
            regulatoryScoreChanged?.endsWith('genie3') &&
            geneSetGenie3Changed
          ) {
            this.fetchAndUpdate(
              regulatoryScoreChanged,
              geneSetGenie3Changed,
            );
          }
        }, 100);
      } else {
        // Clear Genie3 graph if no gene set is selected
        d3.select(graphIdGenie3).selectAll('*').remove();
        compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
      }
    }

    if (geneSetSpongeChanged !== previousGeneSetSpongeChanged) {
      compare ? this.previousGeneSetSpongeCompare = this.selectedGeneSetSponge : this.previousGeneSetSponge = this.selectedGeneSetSponge;
      if (geneSetSpongeChanged) {
        d3.select(graphIdSponge).selectAll('*').remove();
        console.log('Updating Sponge graph for:', geneSetSpongeChanged);
        console.log(
          'Sponge targets available:',
          this.geneSetsSponge[geneSetSpongeChanged]?.length || 0,
        );
        compare ? this.isLoadingSpongeCompare = true : this.isLoadingSponge = true;
        setTimeout(() => {
          this.updateSubgraphSponge(compare);
          if (
            regulatoryScoreChanged?.endsWith('sponge') &&
            geneSetSpongeChanged
          ) {
            this.fetchAndUpdate(
              regulatoryScoreChanged,
              geneSetSpongeChanged,
            );
          }
        }, 100);
      } else {
        // Clear Sponge graph if no gene set is selected
        d3.select(graphIdSponge).selectAll('*').remove();
        compare ? this.isLoadingSpongeCompare = false : this.isLoadingSponge = false;
      }
    }
  }

  public selectGeneSetFromTable(gene: string, networkType: 'genie3' | 'sponge'): void {
    if (networkType === 'genie3') {
      compare ? this.selectedGeneSetGenie3Compare = gene : this.selectedGeneSetGenie3 = gene;
      compare ? this.previousGeneSetGenie3Compare = null : this.previousGeneSetGenie3 = null; // Force update
    } else if (networkType === 'sponge') {
      compare ? this.selectedGeneSetSpongeCompare = gene : this.selectedGeneSetSponge = gene;
      compare ? this.previousGeneSetSpongeCompare = null : this.previousGeneSetSponge = null; // Force update
    }
    // Trigger the graph update
    this.onGeneSetChange(compare);
  }

  public onGeneSelectedFromTable(event: { gene: string; action: string }, networkType: 'genie3' | 'sponge', compare: boolean = false): void {
    // Only select if this action is for the corresponding network type
    if (networkType === 'genie3' && event.action.includes('genie3')) {
      this.selectGeneSetFromTable(event.gene, 'genie3', compare);
    } else if (networkType === 'sponge' && event.action.includes('sponge')) {
      this.selectGeneSetFromTable(event.gene, 'sponge', compare);
    }
  }

  /**
   * Check if initial app loading is complete.
   * Sets isAppInitializing to false once all initial data fetches are done.
   */
  private checkInitializationComplete(compare: boolean = false): void {
    if (compare) {
      if (!this.isLoadingCompare && !this.isLoadingGenie3Compare && !this.isLoadingSpongeCompare && !this.isLoadingRegulatoryScoresCompare) {
        this.isAppInitializing = false;
        console.log('[Init] Compare view initialization complete - hiding loader');
      }
    }
    else {
      if (!this.isLoadingHexagons && !this.isLoadingGenie3 && !this.isLoadingSponge && !this.isLoadingRegulatoryScores) {
        this.isAppInitializing = false;
        console.log('[Init] App initialization complete - hiding loader');
      }
    }
  }

  public selectCellFromCluster(cell: CellFeature, compare: boolean = false): void {
    if (compare) {
      this.selectedCellCompare = cell;
      setTimeout(() => this.updateSubgraphGenie3(true), 0);
    } else {
      this.selectedCell = cell;
      setTimeout(() => this.updateSubgraphGenie3(false), 0);
    }
  }

  public clearClusterData(compare: boolean = false): void {
    this.selectedCluster = null;
    this.clusterCells = [];
    this.clusterCellTypes = [];
    this.coOccurrenceData = []; // Clear co-occurrence data

    if (compare) {
      this.selectedCellCompare = null;
      // If you have compare-specific cluster state, clear it here as well
    } else {
      this.selectedCell = null;
    }

    this.resetClusterExtension();
    this.updateHexColors();
  }

  private calculateClusterStats(compare: boolean = false): void {
    const clusterCells = compare ? this.compareClusterCells : this.clusterCells;
    if (!clusterCells || clusterCells.length === 0) return;

    // Calculate cell type distribution using existing cell_type property
    const cellTypeMap = new Map<string, number>();
    clusterCells.forEach((cell) => {
      const cellType = cell.properties.cell_type;
      cellTypeMap.set(cellType, (cellTypeMap.get(cellType) || 0) + 1);
    });

    const clusterCellTypes = Array.from(cellTypeMap.entries())
      .map(([type, count]) => ({
        type,
        count,
        percentage: ((count / clusterCells.length) * 100).toFixed(1),
      }))
      .sort((a, b) => b.count - a.count);

    if (compare) {
      this.compareClusterCellTypes = clusterCellTypes;
    } else {
      this.clusterCellTypes = clusterCellTypes;
    }

    const clusterAnnotation = this.getLeidenClusterAnnotation(this.selectedCluster!);

    if (clusterCells.length > 0) {
      const centralityAvg = {
        degree_centrality: clusterAnnotation?.centrality?.['degree_centrality'] ?? 0,
        average_clustering: clusterAnnotation?.centrality?.['average_clustering'] ?? 0,
        closeness_centrality: clusterAnnotation?.centrality?.['closeness_centrality'] ?? 0,
      };
      if (compare) {
        this.compareClusterCentralityAvg = centralityAvg;
      } else {
        this.clusterCentralityAvg = centralityAvg;
      }
    }
  }

  private extendCluster(selectedCluster: number, compare: boolean = false): void {
    const features = compare ? this.compareFeatures : this.features;

    this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => {
        // Return original path for non-selected hexagons
        const projection = d3.geoIdentity().fitSize([1200, 1000], {
          type: 'FeatureCollection',
          features: features,
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

  private resetClusterExtension(compare: boolean = false): void {
    const features = compare ? this.compareFeatures : this.features;
    const projection = d3.geoIdentity().fitSize([1200, 1000], {
      type: 'FeatureCollection',
      features: features,
    });
    const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

    targetGroup
      .selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => pathGenerator(d) || '')
      .attr('stroke-width', '1px')
      .attr('stroke', 'transparent')
      .style('opacity', 0.8);

    // Reinitialize the mouseleave event
    targetGroup
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
      this.http.get<{ cell_types: string[] }>(
        `${this.sessionService.apiUrl}/api/cell_types${params}`,
        { withCredentials: true }
      ).subscribe({
        next: resp => { this.availableCellTypes = resp.cell_types ?? []; },
        error: () => { this.availableCellTypes = []; }
      });
    } else {
      this.availableMotifs = [];
      this.availableCellTypes = [];
    }
    console.log('availableCellTypes:', this.availableCellTypes);
  }

  public computeFootprint(): void {
    if (this.footprintMotifs.length === 0) return;
    this.isComputingFootprint = true;
    this.footprintComputeError = '';
    const body = new FormData();
    this.footprintMotifs.forEach(m => body.append('motif', m));
    body.append('cluster_by', this.footprintClusterBy);
    // Pass dataset_id for rescanned-dataset fallback (session may not have adata_path set)
    if (this.selectedDataset?.id) {
      body.append('dataset_id', this.selectedDataset.id);
    }
    this.http.post<{ results: { footprint_url: string; relative_path: string }[] }>(
      `${this.sessionService.apiUrl}/api/compute_footprint`,
      body,
      { withCredentials: true }
    ).subscribe({
      next: resp => {
        this.isComputingFootprint = false;
        const newUrls = (resp.results ?? []).map(r =>
          this.sanitizer.bypassSecurityTrustResourceUrl(
            `${this.sessionService.apiUrl}${r.footprint_url}`
          )
        );
        this.onDemandFootprintUrls = [...this.onDemandFootprintUrls, ...newUrls];
      },
      error: err => {
        this.isComputingFootprint = false;
        this.footprintComputeError =
          err?.error?.detail ?? 'Footprint computation failed. Check the server logs.';
      }
    });
  }

  private renderNhoodHeatmap(compare: boolean = false): void {
    const cell = compare ? this.selectedCellCompare : this.selectedCell;
    const leiden = cell?.properties?.leiden;
    const clusterAnnotation = this.getLeidenClusterAnnotation(leiden, compare);
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
      console.error(`Container ${containerId} not found for rendering heatmap`);
      return;
    }

    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, { displayModeBar: false });
  }

  public closeSidenav(compare: boolean = false): void {
    if (compare) {
      this.selectedCellCompare = null;
    } else {
      this.selectedCell = null;
    }
    this.updateHexColors();
  }

  public updateCoOccurrenceTable(compare: boolean = false): void {
    const selectedCluster = compare ? this.selectedClusterCompare : this.selectedCluster;

    // 1. Clear state immediately
    if (compare) {
      this.compareCoOccurrenceData = [];
    } else {
      this.coOccurrenceData = [];
    }

    if (selectedCluster === null || selectedCluster === undefined) return;

    const clusterAnnotation = this.getLeidenClusterAnnotation(selectedCluster, compare);

    // 2. Validate the 2D array exists for this specific cluster
    if (!clusterAnnotation || !Array.isArray(clusterAnnotation.co_occurrence)) {
      console.warn(`[Co-occurrence] No data matrix found for Cluster ${selectedCluster}`);
      return;
    }

    // This is your [[2.54, 2.10...], [...]] array
    const matrix = clusterAnnotation.co_occurrence;
    const interval = compare ? this.selectedIntervalCompare : this.selectedInterval;
    const data: number[] = [];

    // 3. Safely extract the data
    try {
      // We rely on the actual length of the matrix rather than a hardcoded cluster count
      const numTargetClusters = matrix.length;

      for (let j = 0; j < numTargetClusters; j++) {
        const targetClusterIntervals = matrix[j];

        if (Array.isArray(targetClusterIntervals) && targetClusterIntervals.length > interval) {
          const val = targetClusterIntervals[interval];
          // Guarantee a clean number goes into the array
          data.push(typeof val === 'number' && !isNaN(val) ? val : 0);
        } else {
          data.push(0); // Safe fallback
        }
      }
    } catch (error) {
      console.error('[Co-occurrence] Error parsing 2D matrix:', error);
      data.length = 0; // Wipe on critical failure
    }

    // 4. Update the component state
    if (compare) {
      this.compareCoOccurrenceData = data;
      // Sync the cluster count so the HTML table headers match the data length
      this.clusterCountCompare = data.length;
    } else {
      this.coOccurrenceData = data;
      this.clusterCount = data.length;
    }

    if (data.length > 0) {
      this.calculateCoOccurrenceThreshold(compare);
    }
  }

  private calculateCoOccurrenceThreshold(compare: boolean = false): void {
    const data = compare ? this.compareCoOccurrenceData : this.coOccurrenceData;

    const allValues = data.flat().filter((val) => val > 0);
    if (allValues.length > 0) {
      allValues.sort((a, b) => a - b);
      const percentile75 = Math.floor(allValues.length * 0.75);
      compare ? this.compareCoOccurrenceThreshold = allValues[percentile75] : this.coOccurrenceThreshold = allValues[percentile75];
    }
  }

  public getCoOccurrenceColor(value: number, compare: boolean = false): string {
    if (value === 0) return '#f8f9fa';

    const data = compare ? this.compareCoOccurrenceData : this.coOccurrenceData;
    const selectedCluster = compare ? this.selectedClusterCompare : this.selectedCluster;
    const colorScale = compare ? this.continuousColorScaleCompare : this.continuousColorScale;

    const coOccurrenceWithoutSameCluster = data.filter(
      (_, index) => index !== selectedCluster,
    );
    const maxValue = Math.max(...coOccurrenceWithoutSameCluster);
    const intensity = Math.min(value / maxValue, 1);
    return colorScale(intensity);
  }

  public getIntervalStats(compare: boolean = false): { min: number; max: number; avg: number } {
    const data = compare ? this.compareCoOccurrenceData : this.coOccurrenceData;
    const allValues = data.flat().filter((val) => val > 0);
    if (allValues.length === 0) return { min: 0, max: 0, avg: 0 };

    const min = Math.min(...allValues);
    const max = Math.max(...allValues);
    const avg = allValues.reduce((sum, val) => sum + val, 0) / allValues.length;

    return { min, max, avg: Math.round(avg * 100) / 100 };
  }

  async getRegulatoryScoresforSpots(barcode: string, datasetId?: string, compare: boolean = false) {
    const datasetQuery = datasetId ? `?dataset_id=${encodeURIComponent(datasetId)}` : '';
    this.sessionService.callWithSession(() =>
      this.http.get(
        `${this.sessionService.apiUrl}/obsm/regulatory_scores/cell/${barcode}${datasetQuery}`,
        { withCredentials: true },
      ),
    ).subscribe({
      next: (res) => {
        const rawData = res as { [scoreType: string]: { [element: string]: number } };

        const genie3Data: TableData = {};
        const spongeData: TableData = {};
        const genie3ElementsSet = new Set<string>();
        const spongeElementsSet = new Set<string>();

        for (const [scoreType, scores] of Object.entries(rawData)) {
          if (scoreType.endsWith('_genie3')) {
            genie3Data[scoreType] = scores;
            Object.keys(scores).forEach(element => genie3ElementsSet.add(element));
          } else if (scoreType.endsWith('_sponge')) {
            spongeData[scoreType] = scores;
            Object.keys(scores).forEach(element => spongeElementsSet.add(element));
          }
        }
        if (compare) {
          this.genie3RawDataCompare = genie3Data;
          this.spongeRawDataCompare = spongeData;
          this.genie3ElementsCompare = Array.from(genie3ElementsSet);
          this.spongeElementsCompare = Array.from(spongeElementsSet);
        } else {
          this.genie3RawData = genie3Data;
          this.spongeRawData = spongeData;
          this.genie3Elements = Array.from(genie3ElementsSet);
          this.spongeElements = Array.from(spongeElementsSet);
        }
      },
      error: (err) => {
        if (compare) {
          this.genie3RawDataCompare = {};
          this.spongeRawDataCompare = {};
          this.genie3ElementsCompare = [];
          this.spongeElementsCompare = [];
        } else {
          this.genie3RawData = {};
          this.spongeRawData = {};
          this.genie3Elements = [];
          this.spongeElements = [];
        }
        console.error(
          `[Backend] Failed to load regulatory scores for ${barcode}`,
          err,
        );
      }
    });
  }

  async fetchAndUpdate(columnName: string, index: string, compare: boolean = false) {
    const tokenType = compare ? 'obsm_compare' : 'obsm_main';
    const token = this.nextRequestToken(tokenType);
    const safeIndex = encodeURIComponent(index);
    const isGeneExpression = columnName === 'gene_expression';

    const baseRequest = isGeneExpression
      ? `${this.sessionService.apiUrl}/X/${safeIndex}`
      : `${this.sessionService.apiUrl}/obsm/${encodeURIComponent(columnName)}/${safeIndex}`;

    // Change query dataset based on whether we're updating the main view or compare view

    const datasetQuery = compare
      ? this.selectedDatasetCompare?.id
        ? `?dataset_id=${encodeURIComponent(this.selectedDatasetCompare.id)}`
        : ''
      : this.selectedDataset?.id
        ? `?dataset_id=${encodeURIComponent(this.selectedDataset.id)}`
        : '';

    const request = `${baseRequest}${datasetQuery}`;

    this.sessionService
      .callWithSession(() =>
        this.http.get(
          request,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          if (token !== this.requestTokens[tokenType]) return;
          const data = res as { [barcode: string]: any };

          // 1. Update ONLY the main view if compare is false
          if (!compare && this.features) {
            for (const feature of this.features) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                // HARDCODE the property name. Never use this.selectedView here!
                feature.properties['regulatory_scores'] = data[barcode];
              }
            }
            console.log(`[Backend] Updated main view property 'regulatory_scores' from obsm["${columnName}"][${index}]`);
          }

          // 2. Update ONLY the compare view if compare is true
          if (compare && this.compareMode && this.compareFeatures) {
            for (const feature of this.compareFeatures) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                // HARDCODE the property name. Never use this.selectedCompareView here!
                feature.properties['regulatory_scores'] = data[barcode];
              }
            }
            console.log(`[Backend] Updated compare view property 'regulatory_scores' from obsm["${columnName}"][${index}]`);
          }

          this.repaintBothViews();

          // Mark regulatory scores fetch as complete if during initialization
          if (this.isInitializing) {
            this.isLoadingRegulatoryScores = false;
            this.checkInitializationComplete(compare);
          }
        },
        error: (err) => {
          console.error(
            `[Backend] Failed to load adata.obsm["${columnName}][${index}]`,
            err,
          );
          // Mark as complete even on error
          if (this.isInitializing) {
            this.isLoadingRegulatoryScores = false;
            this.checkInitializationComplete(compare);
          }
        },
      });
  }

  public onRegulatoryScoreChange(compare: boolean = false): void {
    const selectedRegulatoryScore = compare ? this.selectedRegulatoryScoreCompare : this.selectedRegulatoryScore;
    const selectedGeneSetGenie3 = compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3;
    const selectedGeneSetSponge = compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge;

    if (
      selectedRegulatoryScore?.endsWith('genie3') &&
      selectedGeneSetGenie3 &&
      selectedGeneSetSponge
    ) {
      this.fetchAndUpdate(
        selectedRegulatoryScore,
        selectedGeneSetGenie3,
        compare
      );
    } else if (
      selectedRegulatoryScore?.endsWith('sponge') &&
      selectedGeneSetSponge
    ) {
      this.fetchAndUpdate(
        selectedRegulatoryScore,
        selectedGeneSetSponge,
        compare
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
        const value = this.leidenCentralityProps.includes(this.selectedView)
          ? this.getLeidenClusterAnnotation(d.properties.leiden)?.centrality?.[this.selectedView]
          : d.properties?.[this.selectedView];
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
      .on('click', (event, d) => this.displayCellDetails(event, d));
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

interface GeneStatsResponse {
  gene: string;
  dataset_ids: string[];
  global_min: number;
  global_max: number;
  per_dataset: { [datasetId: string]: { min: number; max: number } };
  missing_in: string[];
}
