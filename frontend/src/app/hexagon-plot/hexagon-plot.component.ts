import { Component, OnInit, OnDestroy, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { Data, Router } from '@angular/router';
import { firstValueFrom, map, Observable, Subject, Subscription, takeUntil } from 'rxjs';

// Visualizations
import * as d3 from 'd3';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import Shepherd from 'shepherd.js';

// Services
import { DatasetService } from '../datasets.service';
import { Dataset } from '../datasets.service';
import { SessionService } from '../session.service';
import { GeoDataService } from '../geo-data.service';
import { TranslationService } from '../translation.service';
import { PathsService } from '../paths.service';
import { DEFAULT_PATHS } from '../constants';
import { InfoService } from '../info.service';
import { HexagonViewComponent } from '../hexagon-view/hexagon-view.component';
import { CellInfoPanelComponent } from '../cell-info-panel/cell-info-panel.component';
import { ClusterInfoPanelComponent } from '../cluster-info-panel/cluster-info-panel.component';
import { CoOccurrencePanelComponent } from '../co-occurrence-panel/co-occurrence-panel.component';
import { RegulatoryTablesPanelComponent } from '../regulatory-tables-panel/regulatory-tables-panel.component';
import { ChromvarCorrelationPanelComponent } from '../chromvar-correlation-panel/chromvar-correlation-panel.component';
import { DifferentialMotifActivityPanelComponent } from '../differential-motif-activity-panel/differential-motif-activity-panel.component';
import { FootprintPanelComponent } from '../footprint-panel/footprint-panel.component';
import { DgeaPanelComponent } from '../dgea-panel/dgea-panel.component';
import { RegulatoryScoresPanelComponent } from '../regulatory-scores-panel/regulatory-scores-panel.component';
import { GrnEvaluationPanelComponent } from '../grn-evaluation-panel/grn-evaluation-panel.component';
import { GrnEvaluationOnDemandPanelComponent } from '../grn-evaluation-on-demand-panel/grn-evaluation-on-demand-panel.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';

// Material
import { MatButtonModule } from '@angular/material/button';
import { MatIconModule } from '@angular/material/icon';
import { MatTooltipModule } from '@angular/material/tooltip';
import { MatDialogModule } from '@angular/material/dialog';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatFormField } from '@angular/material/form-field';
import { MatOption } from '@angular/material/autocomplete';
import { MatLabel } from '@angular/material/form-field';
import { MatSelect } from '@angular/material/select';
import { MatExpansionModule } from '@angular/material/expansion';
import { MatTableModule } from '@angular/material/table';
import { MatDividerModule } from '@angular/material/divider';
import { MatTabsModule } from '@angular/material/tabs';
import { MatTabChangeEvent, MatTabGroup, MatTab } from '@angular/material/tabs';
import { DomSanitizer, SafeResourceUrl } from '@angular/platform-browser';
import { MatInputModule } from '@angular/material/input';
import { MatCheckboxModule } from '@angular/material/checkbox';



@Component({
  selector: 'app-hexagon-plot',
  imports: [CommonModule, FormsModule, FilterableTableComponent, HexagonViewComponent, CellInfoPanelComponent, ClusterInfoPanelComponent, CoOccurrencePanelComponent, RegulatoryTablesPanelComponent, ChromvarCorrelationPanelComponent, DifferentialMotifActivityPanelComponent, FootprintPanelComponent, DgeaPanelComponent, RegulatoryScoresPanelComponent, GrnEvaluationPanelComponent, GrnEvaluationOnDemandPanelComponent, MatButtonModule, MatIconModule, MatTooltipModule, MatDialogModule, MatProgressSpinnerModule, MatFormField, MatLabel, MatOption, MatSelect, MatExpansionModule, MatTableModule, MatDividerModule, MatTabsModule, MatInputModule, MatCheckboxModule],
  standalone: true,
  templateUrl: './hexagon-plot.component.html',
  styleUrls: ['./hexagon-plot.component.scss'],
})
export class HexagonPlotComponent implements OnInit, OnDestroy, AfterViewInit {
  @ViewChild('mainTabGroup', { static: false }) tabGroup?: MatTabGroup;
  @ViewChild('compareTabGroup', { static: false }) tabGroupCompare?: MatTabGroup;
  @ViewChild('cellInfoTab', { static: false, read: MatTab }) cellInfoTab?: MatTab;
  @ViewChild('clusterInfoTab', { static: false, read: MatTab }) clusterInfoTab?: MatTab;
  @ViewChild('cellInfoTabCompare', { static: false, read: MatTab }) cellInfoTabCompare?: MatTab;
  @ViewChild('clusterInfoTabCompare', { static: false, read: MatTab }) clusterInfoTabCompare?: MatTab;
  @ViewChild('mainDgeaPanel') mainDgeaPanel?: DgeaPanelComponent;
  @ViewChild('compareDgeaPanel') compareDgeaPanel?: DgeaPanelComponent;
  @ViewChild('mainView') mainView?: HexagonViewComponent;
  @ViewChild('compareView') compareView?: HexagonViewComponent;
  @ViewChild('mainClusterInfo') mainClusterInfo?: ClusterInfoPanelComponent;
  @ViewChild('compareClusterInfo') compareClusterInfo?: ClusterInfoPanelComponent;
  private _resizeHandler: any = null;
  // Removed manual subscription tracker - all subscriptions use takeUntil now
  footprintPlotUrls: SafeResourceUrl[] = [];
  footprintPlotUrlsCompare: SafeResourceUrl[] = [];
  onDemandFootprintUrls: SafeResourceUrl[] = [];
  onDemandFootprintUrlsCompare: SafeResourceUrl[] = [];
  availableMotifs: string[] = [];
  availableMotifsCompare: string[] = [];
  availableCellTypes: string[] = [];
  availableCellTypesCompare: string[] = [];
  footprintMotifs: string[] = [];
  footprintMotifsCompare: string[] = [];
  motifSearchQuery: string = '';
  motifSearchQueryCompare: string = '';
  footprintClusterBy: string = 'cell_type';
  footPrintClusterByCompare: string = 'cell_type';
  isComputingFootprint: boolean = false;
  isComputingFootprintCompare: boolean = false;
  footprintComputeError: string = '';
  footprintComputeErrorCompare: string = '';

  builtinDatasets$: Observable<Dataset[]>;
  uploadedDatasets$: Observable<Dataset[]>;
  tangram_Datasets$: Observable<Dataset[]>;
  selectedDataset: Dataset | null = null;
  selectedDatasetCompare: Dataset | null = null;
  /** Stable bound reference passed to both HexagonViewComponent instances as
   * [propertyAvailableFn], since passing `this.propertyAvailable` directly would lose its
   * `this` binding once invoked from the child. */
  public propertyAvailableBound = (prop: string, isCompare: boolean): boolean =>
    this.propertyAvailable(prop, isCompare);
  /** Same reasoning as propertyAvailableBound, one per side since GrnEvaluationPanelComponent/
   * GrnEvaluationOnDemandPanelComponent's getImageUrl @Input() takes only imageName, not compare. */
  public grnImageUrlBound = (imageName: string): string => this.getGrnImageUrl(imageName, false);
  public grnImageUrlBoundCompare = (imageName: string): string => this.getGrnImageUrl(imageName, true);

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
    public infoService: InfoService,
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

  // svg/g/g_compare/svg_compare/g_paths(*)/baseLayer/detailLayer/currentTransform/
  // currentPathGenerator/detailSize/detailVisible/detailScreenPos/keydownHandler all moved into
  // HexagonViewComponent (accessed here via the mainView/compareView @ViewChild refs where the
  // parent still needs them — see getViewVariablesToUpdate/updateHexColors).

  // ======= Xenium performance state (still parent-owned: computed in loadAndRenderData,
  // read by both the parent's updateHexColors and passed to HexagonViewComponent.renderHexagons) =======
  private fullFeatures: CellFeature[] = [];
  private isXenium = false;


  public selectedCell: CellFeature | null = null;
  public selectedCellCompare: CellFeature | null = null;
  public selectedCluster: number | null = null;
  public selectedClusterCompare: number | null = null;

  public selectedView = 'regulatory_scores';
  public selectedCompareView: string = 'regulatory_scores';

  public selectedItemByView: { [view: string]: string | null } = {};
  public selectedItemByViewCompare: { [view: string]: string | null } = {};

  public regulatoryScoreDisplayMode: 'raw' | 'moranI' | 'gearyC' = 'raw';
  public regulatoryScoreDisplayModeCompare: 'raw' | 'moranI' | 'gearyC' = 'raw';

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

  private geneDomainToken: number = 0;
  private geneDomainCache = new Map<string, { min: number; max: number; expiresAt: number }>();
  private geneDomainCacheTtlMs: number = 10 * 60 * 1000;
  private sharedGeneExpressionDomain: { min: number; max: number } | null = null;
  private sharedGeneExpressionContextKey: string | null = null;
  public dgeaReady: boolean = false;
  public dgeaReadyCompare: boolean = false;

  public selectedInterval: number = 0;
  public selectedIntervalCompare: number = 0;

  public features: CellFeature[] = [];
  public compareFeatures: CellFeature[] = [];

  public meta: { [key: string]: any } = {};
  public metaCompare: { [key: string]: any } = {};

  // Selected groups for the DGEA comparison (bound to the dropdowns)
  public selectedDgeaObsCol: string = 'cell_type';
  public selectedDgeaObsColCompare: string = 'cell_type';
  public selectedDgeaGroup1: string | null = null;
  public selectedDgeaGroup1Compare: string | null = null;
  public selectedDgeaGroup2: string | null = null;
  public selectedDgeaGroup2Compare: string | null = null;
  public dgeaVsAll: boolean = false;
  public dgeaVsAllCompare: boolean = false;


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
  public genie3WidthCompare: number = 600;
  public spongeWidth: number = 600;
  public spongeWidthCompare: number = 600;

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
  private mainMapLoadingCount: number = 0;
  private compareMapLoadingCount: number = 0;
  public isLoadingSponge: boolean = false;
  public isLoadingGenie3: boolean = false;
  public isLoadingGenie3Compare: boolean = false;
  public isLoadingSpongeCompare: boolean = false;

  public isLoadingRegulatoryScores: boolean = false; // Track regulatory scores fetches during init
  public isLoadingRegulatoryScoresCompare: boolean = false; // Track regulatory scores fetches during init for compare view
  public isAppInitializing: boolean = true; // Only for initial app load

  public setMapLoading(loading: boolean, compare: boolean = false): void {
    if (compare) {
      if (loading) {
        this.compareMapLoadingCount++;
      } else {
        this.compareMapLoadingCount = Math.max(0, this.compareMapLoadingCount - 1);
      }
      this.isLoadingCompare = this.compareMapLoadingCount > 0;
    } else {
      if (loading) {
        this.mainMapLoadingCount++;
      } else {
        this.mainMapLoadingCount = Math.max(0, this.mainMapLoadingCount - 1);
      }
      this.isLoadingHexagons = this.mainMapLoadingCount > 0;
    }
  }

  public onTableLoadingChange(event: { loading: boolean; isCompare: boolean }): void {
    this.setMapLoading(event.loading, event.isCompare);
  }

  // Returns true if initial map loading is in progress (ONLY for first time)
  public get isInitializing(): boolean {
    return this.isAppInitializing && this.isLoadingHexagons;
  }

  // Co-occurrence table
  public coOccurrenceData: number[] = [];
  public coOccurrenceThreshold: number = 0.5;
  public maxInterval: number = 49;
  public clusterCount: number = 10;
  // Leiden cluster ids available for this dataset, for the cluster-select dropdown.
  public availableClusterIds: number[] = [];

  // Compare Co-occurrence table
  public compareCoOccurrenceData: number[] = [];
  public compareCoOccurrenceThreshold: number = 0.5;
  public maxIntervalCompare: number = 49;
  public clusterCountCompare: number = 10;
  public availableClusterIdsCompare: number[] = [];

  public currentDataSetSupportsCompare: boolean = true;
  public isXeniumDatasetSelected: boolean = false;

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
  public compareShowGlobalLigandReceptorScores: boolean = true;
  public showMoranI: boolean = true;
  public compareShowMoranI: boolean = true;
  // Categorical scale for cell_type and other non-leiden categorical properties.
  public colorScale = d3.scaleOrdinal<string>(d3.schemeSet2);
  // Dedicated leiden-cluster palette, visually distinct from colorScale above so it's clear
  // at a glance whether the map/legend is showing clusters or another categorical property.
  // Only ever used for the main view's map fill (see HexagonRenderContext.colorScale doc on
  // the pre-existing main/compare cross-contamination this deliberately doesn't touch) — the
  // Cluster Information panel's cell-type swatches keep using colorScale above regardless of
  // the map's current selectedView, since those are always coloring cell_type values.
  public leidenColorScale = d3.scaleOrdinal<string>(d3.schemeTableau10);
  public continuousColorScale = d3.scaleSequential(d3.interpolateYlOrRd).clamp(true);
  // Separate scales for the compare view to avoid cross-contamination
  public colorScaleCompare = d3.scaleOrdinal<string>(d3.schemeSet2);
  public leidenColorScaleCompare = d3.scaleOrdinal<string>(d3.schemeTableau10);
  // Yellow continuous color palette
  public continuousColorScaleCompare = d3.scaleSequential(d3.interpolateYlOrRd).clamp(true);
  public currentLegendDomain: any[] = [];
  public currentLegendType: 'continuous' | 'categorical' = 'categorical';

  // Comparison

  public currentCompareLegendType: 'continuous' | 'categorical' = 'categorical';
  public currentLegendDomainCompare: any[] = [];

  public tfGraphLoading: boolean = false;
  public tfGraphLoadingCompare: boolean = false;
  public selectedGrnView: 'peak_stats' | 'motif_stats' | 'tf_graph' = 'tf_graph';
  public selectedGrnViewCompare: 'peak_stats' | 'motif_stats' | 'tf_graph' = 'tf_graph';

  private tfGraph: tfGraphData | null = null;
  private tfGraphCompare: tfGraphData | null = null;

  // On-demand GRN evaluation
  public onDemandGrnParams = {
    obsKey: 'cell_type',
    cluster: '0',
    name: 'GRN_Evaluation',
    geneSet: '',
    usePriorGrn: false,
  };
  public onDemandGrnParamsCompare = {
    obsKey: 'cell_type',
    cluster: '0',
    name: 'GRN_Evaluation',
    geneSet: '',
    usePriorGrn: false,
  };
  public isComputingOnDemandGrn = false;
  public isComputingOnDemandGrnCompare = false;
  public onDemandGrnError = '';
  public onDemandGrnErrorCompare = '';

  public onDemandOutputPath: string | null = null;
  public onDemandOutputPathCompare: string | null = null;

  // Main View - GRN Evaluation Results
  public onDemandGrnRuns: string[] = [];
  public selectedOnDemandRun: string | null = null;
  public grnGraphFilter: 'prior' | 'extended' | 'full' = 'full';
  public grnGraphData: any = null;
  public grnGraphDataPrior: any = null;
  public grnGraphDataExtended: any = null;
  public grnGraphDataFull: any = null;
  public grnPlots: { [directory: string]: string[] } = {};
  public grnPeakStats: any = null;
  public grnMotifStats: any = null;

  public isLoadingGrnGraph: boolean = false;
  public isLoadingGrnGraphPrior: boolean = false;
  public isLoadingGrnGraphExtended: boolean = false;
  public isLoadingGrnGraphFull: boolean = false;
  public isLoadingGrnRuns: boolean = false;
  public isLoadingGrnPlots: boolean = false;
  public isLoadingGrnStats: boolean = false;
  public grnTabError: string | null = null;

  // Compare View - GRN Evaluation Results (separate properties for independence)
  public onDemandGrnRunsCompare: string[] = [];
  public selectedOnDemandRunCompare: string | null = null;
  public grnGraphFilterCompare: 'prior' | 'extended' | 'full' = 'full';
  public grnGraphDataCompare: any = null;
  public grnGraphDataPriorCompare: any = null;
  public grnGraphDataExtendedCompare: any = null;
  public grnGraphDataFullCompare: any = null;
  public grnPlotsCompare: { [directory: string]: string[] } = {};
  public grnPeakStatsCompare: any = null;
  public grnMotifStatsCompare: any = null;

  public isLoadingGrnRunsCompare: boolean = false;
  public isLoadingGrnGraphCompare: boolean = false;
  public isLoadingGrnGraphPriorCompare: boolean = false;
  public isLoadingGrnGraphExtendedCompare: boolean = false;
  public isLoadingGrnGraphFullCompare: boolean = false;
  public isLoadingGrnPlotsCompare: boolean = false;
  public isLoadingGrnStatsCompare: boolean = false;
  public grnTabErrorCompare: string | null = null;


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
          const dsType = dataset?.dataset_type?.toLowerCase();
          this.currentDataSetSupportsCompare = !(dsType === 'multiome' || dsType === 'xenium');

          // Detect if Xenium dataset is selected
          this.isXeniumDatasetSelected = dsType === 'xenium';

          // Auto-disable compare mode if Xenium is selected
          if (this.isXeniumDatasetSelected && this.compareMode) {
            console.log('Xenium dataset selected - disabling comparison view');
            // Setting compareMode false destroys the *ngIf-gated <app-hexagon-view> compare
            // instance, which takes svg_compare/g_compare (now child-owned) with it — no manual
            // DOM-clearing/field-nulling needed here anymore.
            this.compareMode = false;
            this.datasetService.selectDatasetCompare(null);
          }

          this.dataPath = dataset.geojson_path || '';
          this.features = [];
          this.meta = {};
          this.updatePathsFromDataset(dataset, false);
          // Note: updatePathsFromDataset triggers the pathsService listener,
          // which handles the actual D3 loading for the main view.
        } else {
          this.pathsService.updatePaths({ adataPath: undefined }, false);
          this.dataPath = '';
          this.isXeniumDatasetSelected = false;
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
        if (!paths || (!paths.hexagonPath && !paths.adataPath)) return;

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
            this.mainView?.createHexagonPlot();
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
    // keydownHandler cleanup now happens in HexagonViewComponent's own ngOnDestroy.
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

  // Looked up by id (stable regardless of which component renders the div) rather than
  // @ViewChild, since these containers now live inside RegulatoryScoresPanelComponent's own
  // template — a parent's @ViewChild cannot reach into a child component's template.
  private updateGraphWidths(compare: boolean = false): void {
    const genie3El = document.getElementById(compare ? 'aucell_graph_genie3_compare' : 'aucell_graph_genie3');
    const spongeEl = document.getElementById(compare ? 'aucell_graph_sponge_compare' : 'aucell_graph_sponge');

    if (genie3El) {
      if (compare) this.genie3WidthCompare = genie3El.clientWidth || 600;
      else this.genie3Width = genie3El.clientWidth || 600;
    }
    if (spongeEl) {
      if (compare) this.spongeWidthCompare = spongeEl.clientWidth || 600;
      else this.spongeWidth = spongeEl.clientWidth || 600;
    }
  }

  public onRegulatoryGraphsPanelOpened(compare: boolean = false): void {
    // Wait for lazy expansion-panel content to be attached to the DOM.
    setTimeout(() => {
      this.updateGraphWidths(compare);

      if (compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3) {
        this.updateSubgraphGenie3(compare);
      }

      if (compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge) {
        this.updateSubgraphSponge(compare);
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
      this.mainView?.createHexagonPlot();
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

  public activeTabLabel = 'Cell Information';
  public activeTabLabelCompare = 'Compare - Cell Information';

  public isClusterInfoTabActive(compare: boolean = false): boolean {
    const label = compare ? this.activeTabLabelCompare : this.activeTabLabel;
    return label === 'Cluster Information' || label === 'Compare - Cluster Information';
  }

  public syncClusterHighlight(compare: boolean = false): void {
    const isComp = compare;
    const currentView = isComp ? this.selectedCompareView : this.selectedView;
    const targetView = isComp ? this.compareView : this.mainView;

    if (currentView === 'leiden') {
      let clusterId = isComp ? this.selectedClusterCompare : this.selectedCluster;
      if (clusterId === null) {
        this.autoSelectDefaultCluster(isComp);
        clusterId = isComp ? this.selectedClusterCompare : this.selectedCluster;
      }
      if (clusterId !== null) {
        targetView?.extendCluster(clusterId);
      }
    } else {
      targetView?.resetClusterExtension();
    }
  }

  public getTabLabelForProperty(property: string, compare: boolean = false): string {
    const prefix = compare ? 'Compare - ' : '';
    const prop = (property || '').toLowerCase();

    if (prop === 'leiden') {
      return `${prefix}Cluster Information`;
    }
    if (prop === 'cell_type') {
      return `${prefix}Cell Information`;
    }
    if (prop === 'regulatory_scores' || prop.includes('regulatory_score') || prop.endsWith('_genie3') || prop.endsWith('_sponge')) {
      return `${prefix}Regulatory Scores`;
    }
    if (prop === 'co_occurrence') {
      return `${prefix}Cluster Information`;
    }
    if (prop === 'gene_expression') {
      return `${prefix}Gene Expression`;
    }
    if (prop.includes('ligand_receptor') || prop.startsWith('liana')) {
      return `${prefix}Ligand-Receptor Relationships`;
    }
    if (prop.includes('cell_comp_tf_activity')) {
      return `${prefix}Cell Composition TF Activity`;
    }
    if (prop === 'tf_activity') {
      return `${prefix}TF Activity`;
    }
    if (prop === 'pathway_activity') {
      return `${prefix}Pathway Activity`;
    }
    if (prop === 'dgea') {
      return `${prefix}DGEA`;
    }
    if (prop.includes('chromvar') || prop.includes('moran') || prop.includes('geary')) {
      return compare ? "Compare: ChromVar spatial correlation : Moran's I / Geary's C" : "ChromVar spatial correlation : Moran's I / Geary's C";
    }
    if (prop.includes('diff_motif')) {
      return 'Differential Motif Activity';
    }
    if (prop.includes('footprint')) {
      return 'Footprints';
    }
    if (prop.includes('grn')) {
      return compare ? " Compare - GRN Evaluation" : "GRN Evaluation";
    }
    return `${prefix}Cell Information`;
  }

  public jumpToTabByLabel(tabLabel: string, compare: boolean = false): void {
    const group = compare ? this.tabGroupCompare : this.tabGroup;
    if (!group) return;

    const select = () => {
      const tabs: MatTab[] = (group as any)._tabs?.toArray?.() ?? [];
      const targetTab = tabs.find(t => t.textLabel?.trim().toLowerCase() === tabLabel.trim().toLowerCase());
      if (targetTab) {
        const index = tabs.indexOf(targetTab);
        if (index !== -1) {
          group.selectedIndex = index;
          if (compare) {
            this.activeTabLabelCompare = targetTab.textLabel;
            this.syncClusterHighlight(true);
          } else {
            this.activeTabLabel = targetTab.textLabel;
            this.syncClusterHighlight(false);
          }
        }
      }
    };
    select();
    setTimeout(select, 0);
  }

  /**
   * Selects `tab` within `group`, resolved via the group's own tab list rather than a
   * hardcoded index — the Cell/Cluster Information tabs are added/removed with *ngIf, which
   * shifts every other tab's numeric index whenever they appear or disappear.
   */
  private jumpToTab(group: MatTabGroup | undefined, tab: MatTab | undefined): void {
    if (!group || !tab) return;
    const select = () => {
      const tabs: MatTab[] = (group as any)._tabs?.toArray?.() ?? [];
      const index = tabs.indexOf(tab);
      if (index !== -1) {
        group.selectedIndex = index;
        if (group === this.tabGroup) {
          this.activeTabLabel = tab.textLabel;
          this.syncClusterHighlight(false);
        } else if (group === this.tabGroupCompare) {
          this.activeTabLabelCompare = tab.textLabel;
          this.syncClusterHighlight(true);
        }
      }
    };
    select();
    setTimeout(select, 0);
  }


  private nextRequestToken(graphType: string): number {
    if (!this.requestTokens[graphType]) this.requestTokens[graphType] = 0;
    return ++this.requestTokens[graphType];
  }

  public openInfoPage(fragmentId?: string): void {
    // Simply navigate via router - no full page reload
    // This preserves the component state for when user returns
    this.router.navigate(['/info'], fragmentId ? { fragment: fragmentId } : {});
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

    if (compare) {
      this.activeTabLabelCompare = tabLabel;
    } else {
      this.activeTabLabel = tabLabel;
    }
    this.syncClusterHighlight(compare);

    if (!compare && tabLabel === 'Cluster Information') {
      if (this.selectedView !== 'leiden') {
        this.selectedView = 'leiden';
        this.onColorbyPropertyChange(false);
      }
      let cluster = this.selectedCluster;
      if (cluster === null) {
        this.autoSelectDefaultCluster(false);
        cluster = this.selectedCluster;
      }
      if (cluster !== null) {
        this.mainView?.extendCluster(cluster);
        setTimeout(() => this.mainClusterInfo?.renderNhoodHeatmap(
          cluster,
          this.meta?.['leiden_cluster_annotations'],
        ), 300);
      }
      return;
    }
    if (compare && tabLabel === 'Compare - Cluster Information') {
      if (this.selectedCompareView !== 'leiden') {
        this.selectedCompareView = 'leiden';
        this.onColorbyPropertyChange(true);
      }
      let cluster = this.selectedClusterCompare;
      if (cluster === null) {
        this.autoSelectDefaultCluster(true);
        cluster = this.selectedClusterCompare;
      }
      if (cluster !== null) {
        this.compareView?.extendCluster(cluster);
        setTimeout(() => this.compareClusterInfo?.renderNhoodHeatmap(
          cluster,
          this.metaCompare?.['leiden_cluster_annotations'],
        ), 300);
      }
      return;
    }

    if (!compare && tabLabel === 'DGEA') {
      this.selectedView = 'gene_expression';
      this.dgeaReady = !!this.meta?.['dgea']?.[this.selectedDgeaObsCol];
      if (this.dgeaReady) {
        this.initDgeaSelection(compare);
        setTimeout(() => this.renderDgeaHeatmap(), 100);
      }

    }

    if (compare && tabLabel === 'Compare - DGEA') {
      this.selectedCompareView = 'gene_expression';
      this.dgeaReadyCompare = !!this.metaCompare?.['dgea']?.[this.selectedDgeaObsColCompare];
      if (this.dgeaReadyCompare) {
        this.initDgeaSelection(true);
        setTimeout(() => this.renderDgeaHeatmap(true), 100);
      }

    }

    if (!compare && tabLabel === 'GRN Evaluation') {
      this.selectedGrnView = 'tf_graph';

      // Auto-select first "prob" property or default to cell_type
      const probProperty = Object.keys(this.propertyAvailability).find((prop: string) =>
        prop.toLowerCase().includes('prob') && this.propertyAvailability[prop]
      );
      this.selectedView = probProperty || 'cell_type';

      setTimeout(() => {
        this.onColorbyPropertyChange(false);
        this.loadTfGraph(false);
        this.loadPrecomputedGrnGraph(false);
        this.loadPrecomputedGrnPlots(false);
      }, 300);
      return;
    }

    if (compare && tabLabel.includes('Compare - GRN Evaluation')) {
      this.selectedGrnViewCompare = 'tf_graph';

      // Auto-select first "prob" property or default to cell_type
      const probPropertyCompare = Object.keys(this.propertyAvailabilityCompare).find((prop: string) =>
        prop.toLowerCase().includes('prob') && this.propertyAvailabilityCompare[prop]
      );
      this.selectedCompareView = probPropertyCompare || 'cell_type';

      setTimeout(() => {
        this.onColorbyPropertyChange(true);
        this.loadTfGraph(true);
        this.loadPrecomputedGrnGraph(true);
        this.loadPrecomputedGrnPlots(true);
      }, 100);
      return;
    }

    if (compare && tabLabel === 'GRN Evaluation - On Demand') {
      // Auto-select first "prob" property or default to cell_type
      const probPropertyCompare = Object.keys(this.propertyAvailabilityCompare).find((prop: string) =>
        prop.toLowerCase().includes('prob') && this.propertyAvailabilityCompare[prop]
      );
      this.selectedCompareView = probPropertyCompare || 'cell_type';

      setTimeout(() => {
        this.onColorbyPropertyChange(true);
      }, 100);
      return;
    }



    // Map tab labels to view keys
    const tabMap: { [label: string]: string } = compare
      ? {
        'Compare - Regulatory Scores': 'regulatory_scores',
        'Compare - Gene Expression': 'gene_expression',
        'Compare - Ligand-Receptor Relationships': 'ligand_receptor_relationships',
        'Compare - Cell Composition TF Activity': 'cell_comp_tf_activity_similarity',
        'Compare - TF Activity': 'tf_activity',
        'Compare - Pathway Activity': 'pathway_activity',
        'Compare: ChromVar spatial correlation : Moran\'s I / Geary\'s C': 'chromvar_total_sum',
        'Compare - Differential Motif Activity': 'cell_type',
        'Compare - Footprints': 'cell_type',
        'Compare - DGEA': 'gene_expression'
      }
      : {
        'Regulatory Scores': 'regulatory_scores',
        'Gene Expression': 'gene_expression',
        'TF Activity': 'tf_activity',
        'Pathway Activity': 'pathway_activity',
        'Cell Composition TF Activity': 'cell_comp_tf_activity_similarity',
        'Ligand-Receptor Relationships': 'ligand_receptor_relationships',
        "ChromVar spatial correlation : Moran's I / Geary's C": 'chromvar_total_sum',
        'Differential Motif Activity': 'cell_type',
        'Footprints': 'cell_type',
        'DGEA': 'gene_expression'
      };

    newView = tabMap[tabLabel] || null;

    console.log('[Tab Change]', compare ? '[Compare]' : '[Main]', 'newView=', newView);

    if (newView) {
      this.onTabColorChange(newView, compare);
    }
  }

  /**
   * Handle GRN view selection change
   */
  public onGrnViewChange(view: 'peak_stats' | 'motif_stats' | 'tf_graph', compare: boolean = false): void {
    if (view === 'tf_graph') {
      // Load the TF graph when tf_graph view is selected
      setTimeout(() => this.loadTfGraph(compare), 100);
    }
  }

  /**
   * Check if regulatory scores data is available
   */
  private hasRegulatoryScoresData(isCompare = false): boolean {
    const meta = isCompare ? this.metaCompare : this.meta;
    if (!meta || Object.keys(meta).length === 0) return false;

    const hasGenie3 = meta['global_regulatory_scores_genie3'] &&
      Object.keys(meta['global_regulatory_scores_genie3']).length > 0;
    const hasSponge = meta['global_regulatory_scores_sponge'] &&
      Object.keys(meta['global_regulatory_scores_sponge']).length > 0;

    const hasRegulatoryObsm = (
      isCompare ? this.regulatoryObsmKeysCompare : this.regulatoryObsmKeysMain
    ).length > 0;

    return (hasGenie3 || hasSponge) && hasRegulatoryObsm;
  }

  private hasGeneExpressionData(isCompare = false): boolean {
    const meta = isCompare ? this.metaCompare : this.meta;
    if (!meta || Object.keys(meta).length === 0) return false;

    const hasMoranI = meta['moranI'] && Object.keys(meta['moranI']).length > 0;
    const hasGearyC = meta['gearyC'] && Object.keys(meta['gearyC']).length > 0;

    return (hasMoranI || hasGearyC);
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
  public propertyAvailable(prop: string, isCompare = false): boolean {
    if (prop === 'gene_expression') return this.hasGeneExpressionData(isCompare)
    if (prop === 'regulatory_scores') return this.hasRegulatoryScoresData(isCompare);
    if (prop === 'co_occurrence') return this.hasCoOccurrenceData(isCompare);


    const map = isCompare ? this.propertyAvailabilityCompare : this.propertyAvailability;

    if (map && prop in map) {
      return map[prop];
    }

    return false;
  }

  private hasCoOccurrenceData(isCompare = false): boolean {
    const meta = isCompare ? this.metaCompare : this.meta;

    if (!meta || Object.keys(meta).length === 0 || !meta['leiden_cluster_annotations']) return false;

    const annotations = meta['leiden_cluster_annotations'];

    // Check if at least one cluster has a valid 2D co_occurrence array
    return Object.values(annotations).some((ann: any) =>
      ann &&
      Array.isArray(ann.co_occurrence) &&
      ann.co_occurrence.length > 0 &&
      Array.isArray(ann.co_occurrence[0]) // verify it's a 2D array
    );
  }

  public onCompareMode(): void {
    // Prevent compare mode toggle if Xenium dataset is selected
    if (this.isXeniumDatasetSelected) {
      console.warn('Compare mode not supported for Xenium datasets');
      return;
    }

    this.compareMode = !this.compareMode;

    if (this.compareMode) {
      this.isLoadingCompare = true;
      this.regulatoryObsmKeysCompare = [...this.regulatoryObsmKeysMain];

      // Emit the main dataset to the compare stream.
      // This will automatically trigger the listener in ngOnInit!
      if (this.selectedDataset) {
        this.datasetService.selectDatasetCompare(this.selectedDataset);
      }
      this.refreshSharedGeneExpressionDomain();

    } else {
      // svg_compare/g_compare no longer live here (moved into HexagonViewComponent) — the
      // manual DOM-clearing + field-nulling this used to do is now handled automatically by
      // Angular destroying the *ngIf="compareMode"-gated <app-hexagon-view> instance itself.
      try {
        this.currentLegendDomainCompare = [];
        this.currentCompareLegendType = 'categorical';
        if (this.colorScaleCompare && typeof this.colorScaleCompare.domain === 'function') {
          this.colorScaleCompare.domain([] as any);
        }
        if (this.leidenColorScaleCompare && typeof this.leidenColorScaleCompare.domain === 'function') {
          this.leidenColorScaleCompare.domain([] as any);
        }
        if (this.continuousColorScaleCompare && typeof this.continuousColorScaleCompare.domain === 'function') {
          this.continuousColorScaleCompare.domain([0, 1]);
        }
      } catch (e) { }

      // Clear the service state
      this.compareMapLoadingCount = 0;
      this.isLoadingCompare = false;
      this.datasetService.selectDatasetCompare(null);
      this.refreshSharedGeneExpressionDomain();
    }
  }

  // Selected groups for the DGEA comparison (bound to the dropdowns)
  getSelectedDgeaHeatmap(compare: boolean = false): any | null {
    const cmp = this.getSelectedDgeaComparison(compare);
    if (!cmp || cmp['skipped']) return null;
    return cmp['heatmap_context'] ?? null;
  }

  // Returns all available cell type levels for the dropdown selectors
  getDgeaLevels(compare: boolean = false): string[] {
    const meta = compare ? this.metaCompare : this.meta;
    const obsCol = compare ? this.selectedDgeaObsColCompare : this.selectedDgeaObsCol;
    return meta?.['dgea']?.[obsCol]?.['levels'] ?? [];
  }

  // Returns the map of all DGEA comparisons
  getDgeaComparisonMap(compare: boolean = false): { [key: string]: any } {
    const meta = compare ? this.metaCompare : this.meta;
    const obsCol = compare ? this.selectedDgeaObsColCompare : this.selectedDgeaObsCol;
    return meta?.['dgea']?.[obsCol]?.['comparisons'] ?? {};
  }

  // Re-render the DGEA heatmap when the user changes the group selections
  public onDgeaSelectionChange(compare: boolean = false): void {
    setTimeout(() => this.renderDgeaHeatmap(compare), 0);
  }

  getDgeaObsCols(compare: boolean = false): string[] {
    const meta = compare ? this.metaCompare : this.meta;
    return Object.keys(meta?.['dgea'] ?? {});
  }

  public onDgeaObsColChange(compare: boolean = false): void {
    compare ? this.selectedDgeaGroup1Compare = null : this.selectedDgeaGroup1 = null;
    compare ? this.selectedDgeaGroup2Compare = null : this.selectedDgeaGroup2 = null;
    compare ? this.dgeaVsAllCompare = false : this.dgeaVsAll = false;
    const obsCol = compare ? this.selectedDgeaObsColCompare : this.selectedDgeaObsCol;
    compare
      ? this.dgeaReadyCompare = !!this.metaCompare?.['dgea']?.[obsCol]
      : this.dgeaReady = !!this.meta?.['dgea']?.[obsCol];
    this.initDgeaSelection(compare);
    setTimeout(() => this.renderDgeaHeatmap(compare), 0);
  }

  getDgeaObsColLabel(col: string): string {
    if (col === 'cell_type') return 'Cell type';
    if (col === 'leiden') return 'Leiden';
    return col;
  }

  public hasDgeaData(compare: boolean = false): boolean {
    const meta = compare ? this.metaCompare : this.meta;
    return !!meta?.['dgea'] && Object.keys(meta['dgea']).length > 0;
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
  public shownGeneOnPlotCompare: string | null = null;

  public onDgeaGeneSelected(event: { gene: string; action: string; }, compare: boolean): void {
    if (event.action === 'show_on_plot') {
      this.showDgeaGeneOnMainPlot(event.gene, compare);
    }
  }

  public showDgeaGeneOnMainPlot(gene: string, compare: boolean): void {
    compare ? this.shownGeneOnPlotCompare = gene : this.shownGeneOnPlot = gene;
    compare ? this.selectedCompareView = 'gene_expression' : this.selectedView = 'gene_expression';
    this.onColorbyPropertyChange(compare);
    this.fetchAndUpdate('gene_expression', gene, compare);
  }


  // Delegates to the DgeaPanelComponent instance for this side, imperatively, matching the
  // ClusterInfoPanelComponent/renderNhoodHeatmap pattern (called right after mutating DGEA
  // selection state, in the same tick change detection hasn't refreshed @Input()s yet).
  private renderDgeaHeatmap(compare: boolean = false): void {
    const panel = compare ? this.compareDgeaPanel : this.mainDgeaPanel;
    const obsCol = compare ? this.selectedDgeaObsColCompare : this.selectedDgeaObsCol;
    panel?.renderDgeaHeatmap(this.getSelectedDgeaHeatmap(compare), this.getDgeaObsColLabel(obsCol));
  }

  public onDgeaVsAllChange(compare: boolean = false): void {
    this.initDgeaSelection(compare);
    setTimeout(() => this.renderDgeaHeatmap(compare), 0);
  }

  // Returns all available cell type levels for the dropdown selectors
  getSelectedDgeaComparison(compare: boolean = false): any | null {
    const selectedGroup1 = compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1;
    if (!selectedGroup1) return null;

    const comps = this.getDgeaComparisonMap(compare);

    if (compare ? this.dgeaVsAllCompare : this.dgeaVsAll) {
      const vsAllId = this.makeComparisonId(selectedGroup1, 'all');
      return comps[vsAllId] ?? null;
    }

    const selectedGroup2 = compare ? this.selectedDgeaGroup2Compare : this.selectedDgeaGroup2;
    if (!selectedGroup2) return null;
    if (selectedGroup1 === selectedGroup2) return null;

    const directId = this.makeComparisonId(selectedGroup1, selectedGroup2);
    if (comps[directId]) {
      return comps[directId];
    }

    const reverseId = this.makeComparisonId(selectedGroup2, selectedGroup1);
    const reverseCmp = comps[reverseId];
    if (!reverseCmp) return null;

    return this.flipDgeaComparison(reverseCmp, compare);
  }

  private flipDgeaComparison(cmp: any, compare: boolean = false): any {
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
            compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1,
            compare ? this.selectedDgeaGroup2Compare : this.selectedDgeaGroup2,
            ...cmp.heatmap_context.groups.filter(
              (g: string) => {
                const group1 = compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1;
                const group2 = compare ? this.selectedDgeaGroup2Compare : this.selectedDgeaGroup2;
                return g !== group1 && g !== group2;
              }
            )
          ]
          : cmp.heatmap_context.groups
      }
      : null;

    return {
      ...cmp,
      group1: compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1,
      group2: compare ? this.selectedDgeaGroup2Compare : this.selectedDgeaGroup2,
      n1: cmp.n2,
      n2: cmp.n1,
      name: `${compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1} vs ${compare ? this.selectedDgeaGroup2Compare : this.selectedDgeaGroup2}`,
      table: flippedTable,
      heatmap_context: flippedHeatmap
    };
  }

  // Initialize default selections for the DGEA comparison dropdowns
  initDgeaSelection(compare: boolean): void {
    const levels = this.getDgeaLevels(compare);
    if (!levels.length) return;

    const group1 = compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1;
    if (!group1) {
      compare ? this.selectedDgeaGroup1Compare = levels[0] : this.selectedDgeaGroup1 = levels[0];
    }

    if (compare ? this.dgeaVsAllCompare : this.dgeaVsAll) {
      compare ? this.selectedDgeaGroup2Compare = null : this.selectedDgeaGroup2 = null;
      return;
    }

    const group2 = compare ? this.selectedDgeaGroup2Compare : this.selectedDgeaGroup2;
    if (!group2) {
      const resolvedGroup1 = compare ? this.selectedDgeaGroup1Compare : this.selectedDgeaGroup1;
      const firstDifferent = levels.find(x => x !== resolvedGroup1);
      compare ? this.selectedDgeaGroup2Compare = firstDifferent ?? null : this.selectedDgeaGroup2 = firstDifferent ?? null;
    }
  }


  private loadAndRenderData(dataPath: string, compare: boolean = false): void {
    this.setMapLoading(true, compare);

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
          this.setMapLoading(false, compare);
          this.checkInitializationComplete(compare);
          return;
        }
        console.log('Data loaded from:', fullUrl, data);

        if (!data) {
          throw new Error('Failed to load GeoJSON data');
        }

        // Backend serializes `leiden` as a string (GeoJSON property values round-trip as
        // strings); normalize to number in place so every `properties.leiden === clusterId`
        // comparison downstream (selectCluster, displayClusterDetails, extendCluster, ...)
        // works consistently regardless of whether the caller passes a string or number.
        for (const feature of data.features) {
          const leiden = feature.properties?.leiden;
          if (leiden !== undefined && leiden !== null) {
            feature.properties.leiden = Number(leiden);
          }
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
          this.compareView?.createHexagonPlot();
          this.compareFeatures = data.features;
          // Early preview using whatever pathGenerator/settings a previous render already set up
          // (no-op on the very first-ever compare render, since compareView.currentPathGenerator
          // is only set by renderHexagons — this mirrors the original guard exactly). Superseded
          // moments later by the unified renderHexagons() call below regardless.
          if (this.compareView?.g && this.compareView?.currentPathGenerator) {
            const pathGenerator = this.compareView.currentPathGenerator;
            this.compareView.g
              .style('cursor', 'pointer')
              .style('pointer-events', null)
              .selectAll<SVGPathElement, CellFeature>('path')
              .data(this.compareFeatures, (d: any) => d.properties.barcode)
              .join('path')
              .attr('d', (d: CellFeature) => pathGenerator(d) || '')
              .attr('fill', (d: CellFeature) => {
                const value = this.leidenCentralityProps.includes(this.selectedCompareView)
                  ? this.getLeidenClusterAnnotation(d.properties.leiden, true)?.centrality?.[this.selectedCompareView]
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
              .style('opacity', 0.8);
            // No mouseover/mouseleave/click handlers bound here deliberately — this is a
            // throwaway early preview immediately superseded by the unified renderHexagons()
            // call a few lines below, which (re)binds proper interactivity on the same path
            // elements. Attaching handlers here just to have them overwritten moments later
            // would add complexity for a window too brief for a user to interact with.
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
            if (this.mainView?.svg) {
              if (this.isXenium) {
                (this.mainView.svg as any).on('dblclick.zoom', null);
              }
            }
            // set base
            if (this.isXenium) {
              const target = 10000;
              const step = Math.ceil(this.fullFeatures.length / target);
              this.features = this.fullFeatures.filter((_, i) => i % step === 0);

              // window only for xenium
              if (this.isXenium) {
                this.mainView?.initDetailWindow();
                this.mainView?.hideDetailWindow();
                this.mainView?.bindDetailWindowInteractions();
              }

            } else {
              this.features = this.fullFeatures;
            }
          }



          if (compare) {
            this.metaCompare = data.meta;
            console.log('Compare meta loaded:', this.metaCompare, data.meta);
            const availableObsColsCompare = this.getDgeaObsCols(true);
            if (availableObsColsCompare.length && !availableObsColsCompare.includes(this.selectedDgeaObsColCompare)) {
              this.selectedDgeaObsColCompare = availableObsColsCompare[0];
            }
            this.dgeaReadyCompare = !!this.metaCompare?.['dgea']?.[this.selectedDgeaObsColCompare];
            if (this.dgeaReadyCompare) {
              this.initDgeaSelection(true);
            }

          } else {
            this.meta = data.meta;
            const availableObsCols = this.getDgeaObsCols();
            if (availableObsCols.length && !availableObsCols.includes(this.selectedDgeaObsCol)) {
              this.selectedDgeaObsCol = availableObsCols[0];
            }
            this.dgeaReady = !!this.meta?.['dgea']?.['cell_type'];
            if (this.dgeaReady) {
              this.initDgeaSelection(false);
            }
          }




          const leidenClusterAnnotations = compare ? this.metaCompare?.['leiden_cluster_annotations'] : this.meta?.['leiden_cluster_annotations'];
          if (leidenClusterAnnotations && typeof leidenClusterAnnotations === 'object') {
            const clusterIds = Object.keys(leidenClusterAnnotations).map(Number).sort((a, b) => a - b);
            if (compare) {
              this.clusterCountCompare = clusterIds.length;
              this.availableClusterIdsCompare = clusterIds;
            } else {
              this.clusterCount = clusterIds.length;
              this.availableClusterIds = clusterIds;
            }
          }

          const interval = compare ? this.metaCompare?.['interval'] : this.meta?.['interval'];
          if (Array.isArray(interval) && interval.length > 0) {
            if (compare) {
              this.maxIntervalCompare = interval.length - 1;
            } else {
              this.maxInterval = interval.length - 1;
            }
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
            targetAvailabilityMap[prop] = this.hasRegulatoryScoresData(compare);
          } else if (prop === 'gene_expression') {
            targetAvailabilityMap[prop] = this.hasGeneExpressionData(compare);
          }
          else if (prop === 'co_occurrence') {
            targetAvailabilityMap[prop] = this.hasCoOccurrenceData(compare);

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

        const firstAvailableProperty = compare
          ? this.colorablePropertiesCompare.find((p) => this.propertyAvailable(p, true))
          : this.colorableProperties.find((p) => this.propertyAvailable(p, false));

        if (compare) {
          if (this.propertyAvailable('cell_type', true)) {
            this.selectedCompareView = 'cell_type';
          } else if (this.propertyAvailable('leiden', true)) {
            this.selectedCompareView = 'leiden';
          } else if (this.propertyAvailable('regulatory_scores', true)) {
            this.selectedCompareView = 'regulatory_scores';
          } else if (firstAvailableProperty) {
            this.selectedCompareView = firstAvailableProperty;
          } else {
            // Keep empty when nothing is available.
            this.selectedCompareView = this.colorablePropertiesCompare[0] || '';
          }
          this.currentCompareLegendType = this.isContinuousScale(this.selectedCompareView, this.compareFeatures, true) ? 'continuous' : 'categorical';
        } else {
          if (this.propertyAvailable('cell_type', false)) {
            this.selectedView = 'cell_type';
          } else if (this.propertyAvailable('leiden', false)) {
            this.selectedView = 'leiden';
          } else if (this.propertyAvailable('regulatory_scores', false)) {
            this.selectedView = 'regulatory_scores';
          } else if (firstAvailableProperty) {
            this.selectedView = firstAvailableProperty;
          } else {
            // Keep empty when nothing is available.
            this.selectedView = this.colorableProperties[0] || '';
          }
          this.currentLegendType = this.isContinuousScale(this.selectedView, this.features, false) ? 'continuous' : 'categorical';
        }



        (compare ? this.compareView : this.mainView)?.renderHexagons({
          features: compare ? this.compareFeatures : this.features,
          isXenium: this.isXenium,
          fullFeatures: this.fullFeatures,
          selectedView: compare ? this.selectedCompareView : this.selectedView,
          legendType: compare ? this.currentCompareLegendType : this.currentLegendType,
          // Always the main view's ordinal scale — see HexagonRenderContext.colorScale doc.
          colorScale: this.selectedView === 'leiden' ? this.leidenColorScale : this.colorScale,
          continuousColorScale: compare ? this.continuousColorScaleCompare : this.continuousColorScale,
          leidenClusterAnnotations: (compare ? this.metaCompare : this.meta)?.['leiden_cluster_annotations'],
          leidenCentralityProps: this.leidenCentralityProps,
        });

        this.onColorbyPropertyChange(compare);
        this.autoSelectDefaultCluster(compare);
        setTimeout(() => {
          this.setMapLoading(false, compare);
          this.checkInitializationComplete(compare);
        }, 0);
      })
      .catch((error) => {
        console.error('Error loading or rendering data:', error);
        this.setMapLoading(false, compare);
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

  private extractViewValue(feature: CellFeature, view: string, isCompare = false): unknown {
    if (this.leidenCentralityProps.includes(view)) {
      return this.getLeidenClusterAnnotation(feature.properties.leiden, isCompare)?.centrality?.[view];
    }
    return feature.properties[view];
  }

  private collectFiniteValuesForView(
    features: CellFeature[],
    view: string,
    isCompare = false
  ): number[] {
    return features
      .map((f) => this.toNumber(this.extractViewValue(f, view, isCompare)))
      .filter((n) => Number.isFinite(n));
  }

  private getMinMaxForView(
    features: CellFeature[],
    view: string,
    compare: boolean = false
  ): { min: number; max: number } | null {
    const values = this.collectFiniteValuesForView(features, view, compare);
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

  private getMergedMinMaxForView(view: string): { min: number; max: number } | null {
    // Generalized: always merges min/max for both main and compare for a given view
    const mainFeatures = this.features || [];
    const compareFeatures = this.compareFeatures || [];
    if (!mainFeatures.length && !compareFeatures.length) return null;
    const mainMinMax = mainFeatures.length ? this.getMinMaxForView(mainFeatures, view, false) : null;
    const compareMinMax = compareFeatures.length ? this.getMinMaxForView(compareFeatures, view, true) : null;
    if (!mainMinMax && !compareMinMax) return null;
    const min = [mainMinMax?.min, compareMinMax?.min].filter(x => x !== undefined).reduce((a, b) => Math.min(a as number, b as number));
    const max = [mainMinMax?.max, compareMinMax?.max].filter(x => x !== undefined).reduce((a, b) => Math.max(a as number, b as number));
    return { min, max };
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

    // Fetch pending gene stats if needed
    if (pending.length > 0) {
      this.setMapLoading(true, false);
      if (this.compareMode) this.setMapLoading(true, true);

      try {
        const results = await Promise.allSettled(
          pending.map((r) =>
            firstValueFrom(
              this.sessionService.callWithSession(() =>
                this.http.get<GeneStatsResponse>(
                  `${this.sessionService.apiUrl}/X_stats/${encodeURIComponent(r.gene)}?dataset_ids=${encodeURIComponent(r.datasetId)}`,
                  { withCredentials: true },
                )
              )
            )
          )
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
      } finally {
        this.setMapLoading(false, false);
        if (this.compareMode) this.setMapLoading(false, true);
      }
    }

    // Gather all valid cached domains
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

  public onGeneExpressionSelected(event: { gene: string; action: string }, isCompare: boolean): void {
    if (event.action !== 'gene_expression') return;
    const selectedGene = String(event.gene ?? '').trim();
    if (!selectedGene) return;

    if (!isCompare) {
      this.selectedGeneExpressionMain = selectedGene;
    } else {
      this.selectedGeneExpressionCompare = selectedGene;
    }

    this.refreshSharedGeneExpressionDomain();
  }

  public onItemSelected(event: { gene: string; action: string }, view: string, compare: boolean = false): void {
    const selectedMap = compare ? this.selectedItemByViewCompare : this.selectedItemByView;
    selectedMap[view] = event.gene;
    this.fetchAndUpdate(event.action, event.gene, compare, view);
  }

  public onColorbyPropertyChange(compare: boolean = false): void {
    const colorProp = compare ? this.selectedCompareView : this.selectedView;
    const targetView = compare ? this.compareView : this.mainView;
    targetView?.setCurrentView(colorProp);

    // Automatically jump sidebar to matching tab for the selected property
    const matchingTabLabel = this.getTabLabelForProperty(colorProp, compare);
    this.jumpToTabByLabel(matchingTabLabel, compare);

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

    if (colorProp === 'leiden') {
      const cluster = compare ? this.selectedClusterCompare : this.selectedCluster;
      if (cluster !== null) {
        targetView?.extendCluster(cluster);
      }
    } else {
      targetView?.resetClusterExtension();
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
        const clusterAnnotation = this.getLeidenClusterAnnotation(f.properties.leiden, compare);
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

    return allNumbers && !shouldTreatAsCategorical && numericValues.length > 0;
  }


  private getViewVariablesToUpdate(containerName: string) {
    const isMainView = containerName === '#hexbin';
    const viewToUpdate = isMainView ? this.selectedView : this.selectedCompareView;
    const featuresToUpdate = isMainView ? this.features : (this.compareFeatures || []);
    const gToUpdate = isMainView ? this.mainView?.g : this.compareView?.g;
    const ordinalScaleToUpdate = viewToUpdate === 'leiden'
      ? (isMainView ? this.leidenColorScale : this.leidenColorScaleCompare)
      : (isMainView ? this.colorScale : this.colorScaleCompare);
    const continuousScaleToUpdate = isMainView ? this.continuousColorScale : this.continuousColorScaleCompare;
    // Pass the compare flag to isContinuousScale for correct logic
    const isContinuous = this.isContinuousScale(viewToUpdate, featuresToUpdate, !isMainView);
    const legendContainerName = isMainView ? 'svg-legend' : 'svg-legend-compare';
    const legendGradientName = isMainView ? 'svg-legend-gradient' : 'svg-legend-gradient-compare';
    return {
      isMainView,
      view: viewToUpdate,
      svg: isMainView ? this.mainView?.svg : this.compareView?.svg,
      features: featuresToUpdate,
      g: gToUpdate,
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

      if (!hasProperty) {
        console.warn('[updateHexColors] Property not found in features! Available:', Object.keys(firstFeature.properties || {}));
        // Optionally, show a user-facing message here
      }
    }

    const targetView = isMainView ? this.mainView : this.compareView;
    targetView?.setCurrentView(viewToUse);

    const layerToColor = (
      (containerName === '#hexbin' && this.isXenium)
        ? this.mainView?.baseLayer
        : viewVariablesToUpdate.g
    ) as unknown as d3.Selection<SVGGElement, any, any, any>;

    // Safety check: if layerToColor is null/undefined, skip the color update
    if (!layerToColor) {
      console.warn(`[updateHexColors] Layer not found for ${containerName}`);
      return;
    }

    const sel = layerToColor
      .selectAll('path')
      .data(features);

    if (containerName === '#hexbin-compare' && this.compareMode && sel.size() === 0) {
      console.warn('[updateHexColors] Compare layer has no paths; re-initializing compare plot');
      this.loadAndRenderData(this.selectedDatasetCompare?.geojson_path || '', true);
      return;
    }

    valuesRaw = features.map((f) => this.extractViewValue(f, viewToUse, !isMainView));
    const numericValues = valuesRaw.map((v) => this.toNumber(v));

    // Determine if the property is continuous
    const isContinuous = this.isContinuousScale(viewToUse, features, !isMainView);

    if (isContinuous) {
      const compareSharedDomain = this.getPairedContinuousDomainForCompare();
      console.log('[updateHexColors] Computed compareSharedDomain:', compareSharedDomain);
      const sharedDomain = compareSharedDomain ?? (viewToUse === 'gene_expression'
        ? this.getSharedDomainForGeneExpressionView()
        : null);

      let min = sharedDomain ? sharedDomain.min : Math.min(...numericValues);
      let max = sharedDomain ? sharedDomain.max : Math.max(...numericValues);

      if (!Number.isFinite(min) || !Number.isFinite(max)) {
        min = 0;
        max = 1;
      }

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
        .interrupt()
        .attr('fill', (d) => {
          const raw = this.extractViewValue(d, viewToUse, !isMainView);
          const n = this.toNumber(raw);
          return Number.isFinite(n)
            ? viewVariablesToUpdate.continuous((n - min) / (max - min))
            : '#ccc';
        });
    } else {
      // categorical scale
      const domain = [...new Set(valuesRaw.map((v: any) => String(v)))];

      viewVariablesToUpdate.ordinal.domain(domain);
      viewVariablesToUpdate.setLegendDomain(domain);
      viewVariablesToUpdate.setLegendType('categorical');

      sel
        .interrupt()
        .attr('fill', (d) => {
          const raw = this.extractViewValue(d, viewToUse, !isMainView);
          return viewVariablesToUpdate.ordinal(String(raw));
        });
    }

    this.syncClusterHighlight(!isMainView);
    targetView?.updateSelectionHighlight();

    if (containerName === '#hexbin') {
      this.mainView?.refreshDetailWindowIfVisible();
    }
    this.renderLegend(containerName);

  }

  public updateSubgraphGenie3(compare: boolean = false): void {
    const token = this.nextRequestToken('genie3');
    console.log('Updating AUCELL graph for Genie3...');
    compare ? this.isLoadingGenie3Compare = true : this.isLoadingGenie3 = true;
    const geneSet = compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3;
    const datasetId = compare ? this.selectedDatasetCompare?.id : this.selectedDataset?.id;
    const graphContainerId = compare ? '#aucell_graph_genie3_compare' : '#aucell_graph_genie3';
    d3.select(graphContainerId).selectAll('*').remove();

    if (!geneSet || !datasetId) {
      compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
      this.checkInitializationComplete(compare);
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
          if (token !== this.requestTokens['genie3']) {
            compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
            this.checkInitializationComplete(compare);
            return;
          }
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
          if (token !== this.requestTokens['genie3']) {
            compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
            this.checkInitializationComplete(compare);
            return;
          }
          console.error(
            `[Backend] Failed to load Genie3 Connections for["${geneSet}"]`,
            err,
          );
          compare ? this.isLoadingGenie3Compare = false : this.isLoadingGenie3 = false;
          this.checkInitializationComplete(compare);
        },
      });
  }

  private getPairedContinuousDomainForCompare(): { min: number; max: number } | null {
    if (!this.compareMode) {
      console.log('[domain:getPaired] skipped: compare mode is off');
      return null;
    }

    const mainView = this.selectedView;
    const compareView = this.selectedCompareView;

    if (mainView !== 'gene_expression' || compareView !== 'gene_expression') {
      return null;
    }

    const mainFeatures = this.features || [];
    const compareFeatures = this.compareFeatures || [];

    if (!mainFeatures.length || !compareFeatures.length) {
      return null;
    }

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
    const compareMinMax = this.getMinMaxForView(compareFeatures, compareView, true);
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
    slicedEdges = candidateEdges.slice(0, compare ? this.genie3MinEdgesCompare : this.genie3MinEdges);

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
    const width = compare ? this.genie3WidthCompare : this.genie3Width;
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

    setTimeout(() => this.updateGraphWidths(compare), 50);

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
      this.checkInitializationComplete(compare);
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
          if (token !== this.requestTokens['sponge']) {
            compare ? this.isLoadingSpongeCompare = false : this.isLoadingSponge = false;
            this.checkInitializationComplete(compare);
            return;
          }
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
              if (compare) {
                this.spongeSliderDataCompare = {
                  step: sliderData.step || 0.01,
                  min_border: sliderData.min_border || 0,
                  max_border: sliderData.max_border || 100,
                  default_value: sliderData.default_value || 50,
                };
                this.spongePValueCutoffCompare = this.spongeSliderDataCompare.default_value;
              } else {
                this.spongeSliderData = {
                  step: sliderData.step || 0.01,
                  min_border: sliderData.min_border || 0,
                  max_border: sliderData.max_border || 100,
                  default_value: sliderData.default_value || 50,
                };
                this.spongePValueCutoff = this.spongeSliderData.default_value;
              }
            }

            console.log('Sponge Network:', data);

            compare ? this.spongeNetworkCompare = data.map((d) => ({
              source: d.geneA,
              target: d.geneB,
              p_adjusted: d['p.adj'],
              mscore: d['mscor'],
            })) : this.spongeNetwork = data.map((d) => ({
              source: d.geneA,
              target: d.geneB,
              p_adjusted: d['p.adj'],
              mscore: d['mscor'],
            }));

            this.visualizeSpongeSubgraph(compare);
          }
        },
        error: (err) => {
          if (token !== this.requestTokens['sponge']) {
            compare ? this.isLoadingSpongeCompare = false : this.isLoadingSponge = false;
            this.checkInitializationComplete(compare);
            return;
          }
          console.error(
            `[Backend] Failed to load Sponge Connections for["${geneSet}"]`,
            err,
          );
          compare ? this.isLoadingSpongeCompare = false : this.isLoadingSponge = false;
          this.checkInitializationComplete(compare);
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
    let targets = this.geneSetsSponge[regulator] || [];


    let nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    let edges: { source: string; target: string; p_adjusted: number }[] = [];

    let candidateEdges: { source: string; target: string; p_adjusted: number }[] = [];
    let slicedEdges: { source: string; target: string; p_adjusted: number }[] = [];

    if (compare) {
      candidateEdges = this.spongeNetworkCompare.filter((edge) => edge.p_adjusted < this.spongePValueCutoffCompare).map((e) => ({
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
    slicedEdges = candidateEdges.slice(0, compare ? this.spongeMinEdgesCompare : this.spongeMinEdges);

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
    const width = compare ? this.spongeWidthCompare : this.spongeWidth;
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

    setTimeout(() => this.updateGraphWidths(compare), 50);

    console.log('Network visualization complete');
  }


  public onSliderChangeSponge(compare: boolean = false): void {
    if (compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge) {
      compare ? this.isLoadingSpongeCompare = true : this.isLoadingSponge = true;
      setTimeout(() => this.visualizeSpongeSubgraph(compare), 50);
    }
  }


  public onSliderChangeGenie3(compare: boolean = false): void {
    if (compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3) {
      compare ? this.isLoadingGenie3Compare = true : this.isLoadingGenie3 = true;
      setTimeout(() => this.visualizeGenie3Subgraph(compare), 50);
    }
  }

  public analyzeGeneSetInGProfiler(networkType: 'genie3' | 'sponge' = 'genie3'): void {
    const geneSets = networkType === 'sponge' ? this.geneSetsSponge : this.geneSetsGenie3;
    const selectedGeneSet = networkType === 'sponge' ? this.selectedGeneSetSponge : this.selectedGeneSetGenie3;

    if (!selectedGeneSet || !geneSets) {
      console.warn(`No ${networkType} gene set selected`);
      return;
    }

    const regulator = selectedGeneSet;
    const targets = geneSets[regulator] || [];

    const allGenes = [regulator, ...targets];

    console.log(`Analyzing all ${networkType} genes in gProfiler:`, allGenes);

    const gProfilerUrl = this.generateGProfilerUrl(allGenes);

    if (gProfilerUrl) {
      // Open in new tab/window
      window.open(gProfilerUrl, '_blank');
    } else {
      console.warn(
        `Could not generate gProfiler URL for ${networkType} gene set:`,
        selectedGeneSet,
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


  public displayCellDetails(event: MouseEvent, cell: CellFeature, compare: boolean = false): void {
    if (!compare) {
      // 1. Always select the clicked cell
      this.selectedCell = cell;

      // 2. Always update cluster information in the background if leiden is present
      const rawLeiden = cell.properties?.leiden;
      const clusterId = rawLeiden !== undefined && rawLeiden !== null && !isNaN(Number(rawLeiden))
        ? Number(rawLeiden)
        : null;

      if (clusterId !== null) {
        this.selectedCluster = clusterId;
        this.clusterCells = (this.features || []).filter(
          (c) => Number(c.properties?.leiden) === clusterId,
        );
        this.calculateClusterStats(false);
        this.updateCoOccurrenceTable(false);
        setTimeout(() => this.updateSubgraphGenie3(false), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
        setTimeout(() => this.mainClusterInfo?.renderNhoodHeatmap(
          clusterId,
          this.meta?.['leiden_cluster_annotations'],
        ), 300);
      }

      // 3. Fetch regulatory scores if in that view
      if (this.selectedView === 'regulatory_scores') {
        this.getRegulatoryScoresforSpots(cell.properties.barcode, this.selectedDataset?.id, false);
      }

      // 4. View-dependent highlighting & tab navigation
      if (this.selectedView === 'leiden' && clusterId !== null) {
        // Outline entire cluster (thin) + clicked cell (thick)
        this.mainView?.extendCluster(clusterId);
        this.jumpToTab(this.tabGroup, this.clusterInfoTab);
      } else {
        // Outline only clicked cell (thick)
        this.mainView?.resetClusterExtension();
        this.jumpToTab(this.tabGroup, this.cellInfoTab);
      }
    } else {
      // 1. Always select the clicked cell in compare view
      this.selectedCellCompare = cell;

      // 2. Always update compare cluster information in the background if leiden is present
      const rawLeiden = cell.properties?.leiden;
      const clusterId = rawLeiden !== undefined && rawLeiden !== null && !isNaN(Number(rawLeiden))
        ? Number(rawLeiden)
        : null;

      if (clusterId !== null) {
        this.selectedClusterCompare = clusterId;
        this.compareClusterCells = (this.compareFeatures || []).filter(
          (c) => Number(c.properties?.leiden) === clusterId,
        );
        this.calculateClusterStats(true);
        this.updateCoOccurrenceTable(true);
        setTimeout(() => this.updateSubgraphGenie3(true), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDatasetCompare, true), 100);
        setTimeout(() => this.compareClusterInfo?.renderNhoodHeatmap(
          clusterId,
          this.metaCompare?.['leiden_cluster_annotations'],
        ), 300);
      }

      // 3. Fetch regulatory scores if in that view
      if (this.selectedCompareView === 'regulatory_scores') {
        this.getRegulatoryScoresforSpots(
          cell.properties.barcode,
          this.selectedDatasetCompare?.id ?? this.selectedDataset?.id,
          true,
        );
      }

      // 4. View-dependent highlighting & tab navigation
      if (this.selectedCompareView === 'leiden' && clusterId !== null) {
        // Outline entire cluster (thin) + clicked cell (thick)
        this.compareView?.extendCluster(clusterId);
        this.jumpToTab(this.tabGroupCompare, this.clusterInfoTabCompare);
      } else {
        // Outline only clicked cell (thick)
        this.compareView?.resetClusterExtension();
        this.jumpToTab(this.tabGroupCompare, this.cellInfoTabCompare);
      }
    }
  }

  /**
   * Auto-selects the first cluster that has co-occurrence data in the background,
   * without applying map selection borders or selecting an individual cell.
   */
  private autoSelectDefaultCluster(compare: boolean = false): void {
    const meta = compare ? this.metaCompare : this.meta;
    const annotations = meta?.['leiden_cluster_annotations'] as Record<string, any> | undefined;
    if (!annotations) return;

    const defaultClusterId = Object.keys(annotations).find((id) => {
      const co_occurrence = annotations[id]?.co_occurrence;
      return Array.isArray(co_occurrence) && co_occurrence.length > 0 && Array.isArray(co_occurrence[0]);
    });

    if (defaultClusterId !== undefined) {
      const clusterNum = Number(defaultClusterId);
      if (compare) {
        this.selectedClusterCompare = clusterNum;
        this.compareClusterCells = (this.compareFeatures || []).filter(
          (c) => Number(c.properties?.leiden) === clusterNum,
        );
        this.calculateClusterStats(true);
        this.updateCoOccurrenceTable(true);
        setTimeout(() => this.updateSubgraphGenie3(true), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDatasetCompare, true), 100);
      } else {
        this.selectedCluster = clusterNum;
        this.clusterCells = (this.features || []).filter(
          (c) => Number(c.properties?.leiden) === clusterNum,
        );
        this.calculateClusterStats(false);
        this.updateCoOccurrenceTable(false);
        setTimeout(() => this.updateSubgraphGenie3(false), 100);
        setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
      }
    }
  }

  public displayClusterDetails(clusterId: number): void {
    this.selectCluster(clusterId, false);
  }

  public selectCluster(clusterId: number, compare: boolean = false): void {
    if (compare) {
      // 1. Switch colored-by view to Leiden if not already
      if (this.selectedCompareView !== 'leiden') {
        this.selectedCompareView = 'leiden';
        this.onColorbyPropertyChange(true);
      }

      // 2. Set active cluster state
      this.selectedClusterCompare = clusterId;
      this.compareClusterCells = (this.compareFeatures || []).filter(
        (cell) => Number(cell.properties?.leiden) === clusterId,
      );
      this.calculateClusterStats(true);
      this.updateCoOccurrenceTable(true);

      // 3. Clear single cell selection if it was not in this cluster (no auto-selected cell)
      if (this.selectedCellCompare && Number(this.selectedCellCompare.properties?.leiden) !== clusterId) {
        this.selectedCellCompare = null;
      }

      // 4. Highlight the cluster with thin border on the map
      this.compareView?.extendCluster(clusterId);

      setTimeout(() => this.updateSubgraphGenie3(true), 100);
      setTimeout(() => this.renderFootprintPlots(this.selectedDatasetCompare, true), 100);
      setTimeout(() => this.compareClusterInfo?.renderNhoodHeatmap(
        clusterId,
        this.metaCompare?.['leiden_cluster_annotations'],
      ), 300);
    } else {
      // 1. Switch colored-by view to Leiden if not already
      if (this.selectedView !== 'leiden') {
        this.selectedView = 'leiden';
        this.onColorbyPropertyChange(false);
      }

      // 2. Set active cluster state
      this.selectedCluster = clusterId;
      this.clusterCells = (this.features || []).filter(
        (cell) => Number(cell.properties?.leiden) === clusterId,
      );
      this.calculateClusterStats(false);
      this.updateCoOccurrenceTable(false);

      // 3. Clear single cell selection if it was not in this cluster (no auto-selected cell)
      if (this.selectedCell && Number(this.selectedCell.properties?.leiden) !== clusterId) {
        this.selectedCell = null;
      }

      // 4. Highlight the cluster with thin border on the map
      this.mainView?.extendCluster(clusterId);

      setTimeout(() => this.updateSubgraphGenie3(false), 100);
      setTimeout(() => this.renderFootprintPlots(this.selectedDataset), 100);
      setTimeout(() => this.mainClusterInfo?.renderNhoodHeatmap(
        clusterId,
        this.meta?.['leiden_cluster_annotations'],
      ), 300);
    }
  }

  public onGeneSetChange(compare: boolean = false): void {
    const geneSetGenie3Changed = compare ? this.selectedGeneSetGenie3Compare : this.selectedGeneSetGenie3;
    const geneSetSpongeChanged = compare ? this.selectedGeneSetSpongeCompare : this.selectedGeneSetSponge;
    const previousGeneSetGenie3Changed = compare ? this.previousGeneSetGenie3Compare : this.previousGeneSetGenie3;
    const previousGeneSetSpongeChanged = compare ? this.previousGeneSetSpongeCompare : this.previousGeneSetSponge;
    const regulatoryScoreChanged = compare ? this.selectedRegulatoryScoreCompare : this.selectedRegulatoryScore;
    const graphIdGenie3 = compare ? '#aucell_graph_genie3_compare' : '#aucell_graph_genie3';
    const graphIdSponge = compare ? '#aucell_graph_sponge_compare' : '#aucell_graph_sponge';

    if (geneSetGenie3Changed !== previousGeneSetGenie3Changed) {
      compare ? this.previousGeneSetGenie3Compare = this.selectedGeneSetGenie3Compare : this.previousGeneSetGenie3 = this.selectedGeneSetGenie3;
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
              compare
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
      compare ? this.previousGeneSetSpongeCompare = this.selectedGeneSetSpongeCompare : this.previousGeneSetSponge = this.selectedGeneSetSponge;
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
              compare
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

  public selectGeneSetFromTable(gene: string, networkType: 'genie3' | 'sponge', compare: boolean = false): void {
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
    if (networkType === 'genie3' && event.action.includes('genie3')) {
      if (compare) {
        this.selectedRegulatoryScoreCompare = event.action;
      } else {
        this.selectedRegulatoryScore = event.action;
      }

      this.selectGeneSetFromTable(event.gene, 'genie3', compare);
    } else if (networkType === 'sponge' && event.action.includes('sponge')) {

      if (compare) {
        this.selectedRegulatoryScoreCompare = event.action;
      } else {
        this.selectedRegulatoryScore = event.action;
      }
      this.selectGeneSetFromTable(event.gene, 'sponge', compare);
    }
  }

  /**
   * Check if initial app loading is complete.
   * Sets isAppInitializing to false once map data fetch and render are done.
   */
  private checkInitializationComplete(compare: boolean = false): void {
    if (compare) {
      if (!this.isLoadingCompare) {
        this.isAppInitializing = false;
        console.log('[Init] Compare view initialization complete - hiding loader');
      }
    }
    else {
      if (!this.isLoadingHexagons) {
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
    if (compare) {
      this.selectedClusterCompare = null;
      this.compareClusterCells = [];
      this.clusterCellTypesCompare = [];
      this.compareCoOccurrenceData = [];
      this.selectedCellCompare = null;
    } else {
      this.selectedCluster = null;
      this.clusterCells = [];
      this.clusterCellTypes = [];
      this.coOccurrenceData = [];
      this.selectedCell = null;
    }

    (compare ? this.compareView : this.mainView)?.resetClusterExtension(compare ? this.compareFeatures : this.features);
    this.updateHexColors(compare ? '#hexbin-compare' : '#hexbin');
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

    const clusterAnnotation = this.getLeidenClusterAnnotation(
      compare ? this.selectedClusterCompare! : this.selectedCluster!,
      compare
    );

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


  private renderFootprintPlots(dataset: Dataset | null, compare: boolean = false): void {
    const raw = dataset?.footprint_pdf_paths;

    const precomputedPaths =
      Array.isArray(raw)
        ? raw
        : raw && typeof raw === 'object'
          ? Object.values(raw)
          : (dataset?.footprint_list ?? []);

    const urls = precomputedPaths.map(path =>
      this.sanitizer.bypassSecurityTrustResourceUrl(
        `${this.sessionService.apiUrl}/api/download/${path}`
      )
    );
    if (compare) {
      this.footprintPlotUrlsCompare = urls;
    } else {
      this.footprintPlotUrls = urls;
    }

    if (dataset) {
      const params = dataset.id ? `?dataset_id=${encodeURIComponent(dataset.id)}` : '';
      this.http.get<{ motifs: string[] }>(
        `${this.sessionService.apiUrl}/api/motifs${params}`,
        { withCredentials: true }
      ).subscribe({
        next: resp => {
          compare
            ? this.availableMotifsCompare = resp.motifs ?? []
            : this.availableMotifs = resp.motifs ?? [];
        },
        error: () => { compare ? this.availableMotifsCompare = [] : this.availableMotifs = []; }
      });

      this.http.get<{ cell_types: string[] }>(
        `${this.sessionService.apiUrl}/api/cell_types${params}`,
        { withCredentials: true }
      ).subscribe({
        next: resp => {
          compare
            ? this.availableCellTypesCompare = resp.cell_types ?? []
            : this.availableCellTypes = resp.cell_types ?? [];
        },
        error: () => { compare ? this.availableCellTypesCompare = [] : this.availableCellTypes = []; }
      });
    } else if (compare) {
      this.availableMotifsCompare = [];
      this.availableCellTypesCompare = [];
    } else {
      this.availableMotifs = [];
      this.availableCellTypes = [];
    }
  }

  public computeFootprint(compare: boolean = false): void {
    if (compare) {
      if (this.footprintMotifsCompare.length === 0) return;
      this.isComputingFootprintCompare = true;
      this.footprintComputeErrorCompare = '';
      const bodyCompare = new FormData();
      this.footprintMotifsCompare.forEach(m => bodyCompare.append('motif', m));
      bodyCompare.append('cluster_by', this.footPrintClusterByCompare);
      if (this.selectedDatasetCompare?.id) {
        bodyCompare.append('dataset_id', this.selectedDatasetCompare.id);
      }
      this.http.post<{ results: { footprint_url: string; relative_path: string }[] }>(
        `${this.sessionService.apiUrl}/api/compute_footprint`,
        bodyCompare,
        { withCredentials: true }
      ).subscribe({
        next: resp => {
          this.isComputingFootprintCompare = false;
          const newUrls = (resp.results ?? []).map(r =>
            this.sanitizer.bypassSecurityTrustResourceUrl(
              `${this.sessionService.apiUrl}${r.footprint_url}`
            )
          );
          this.onDemandFootprintUrlsCompare = [...this.onDemandFootprintUrlsCompare, ...newUrls];
        },
        error: err => {
          this.isComputingFootprintCompare = false;
          this.footprintComputeErrorCompare =
            err?.error?.detail ?? 'Footprint computation failed. Check the server logs.';
        }
      });
      return;
    }
    else {
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

  public getIntervalBoundaries(compare: boolean = false): number[] {
    const raw = compare ? this.metaCompare?.['interval'] : this.meta?.['interval'];
    return Array.isArray(raw)
      ? raw.map((v: any) => Number(v)).filter((v: number) => Number.isFinite(v))
      : [];
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

  async fetchAndUpdate(columnName: string, index: string, compare: boolean = false, view?: string) {
    const tokenType = compare ? 'obsm_compare' : 'obsm_main';
    const token = this.nextRequestToken(tokenType);
    const safeIndex = encodeURIComponent(index);
    const isGeneExpression = columnName === 'gene_expression';

    // Determine the property name to update based on the view
    // If view is provided, use it; otherwise default to 'regulatory_scores' for backward compatibility
    const propertyToUpdate =
      view || (isGeneExpression ? 'gene_expression' : 'regulatory_scores');

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

    this.setMapLoading(true, compare);

    this.sessionService
      .callWithSession(() =>
        this.http.get(
          request,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          if (token !== this.requestTokens[tokenType]) {
            this.setMapLoading(false, compare);
            return;
          }
          const data = res as { [barcode: string]: any };

          // 1. Update ONLY the main view if compare is false
          if (!compare && this.features) {
            for (const feature of this.features) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                feature.properties[propertyToUpdate] = data[barcode];
              }
            }
            console.log(`[Backend] Updated main view property '${propertyToUpdate}' from obsm["${columnName}"][${index}]`);
          }

          // 2. Update ONLY the compare view if compare is true
          if (compare && this.compareMode && this.compareFeatures) {
            for (const feature of this.compareFeatures) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                feature.properties[propertyToUpdate] = data[barcode];
              }
            }
            console.log(`[Backend] Updated compare view property '${propertyToUpdate}' from obsm["${columnName}"][${index}]`);
          }

          this.repaintBothViews();
          this.setMapLoading(false, compare);

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
          this.setMapLoading(false, compare);
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

    if (selectedRegulatoryScore?.endsWith('genie3') && selectedGeneSetGenie3) {
      this.fetchAndUpdate(selectedRegulatoryScore, selectedGeneSetGenie3, compare);
    } else if (selectedRegulatoryScore?.endsWith('sponge') && selectedGeneSetSponge) {
      this.fetchAndUpdate(selectedRegulatoryScore, selectedGeneSetSponge, compare);
    }
  }


  public onRegulatoryDisplayModeChange(compare: boolean = false): void {
    const displayMode = compare ? this.regulatoryScoreDisplayModeCompare : this.regulatoryScoreDisplayMode;
    const meta = compare ? this.metaCompare : this.meta;

    // Validate that data exists for the selected display mode
    const genie3Key = displayMode === 'raw'
      ? 'global_regulatory_scores_genie3'
      : displayMode === 'moranI'
        ? 'global_regulatory_moranI_genie3'
        : 'global_regulatory_gearyC_genie3';

    const spongeKey = displayMode === 'raw'
      ? 'global_regulatory_scores_sponge'
      : displayMode === 'moranI'
        ? 'global_regulatory_moranI_sponge'
        : 'global_regulatory_gearyC_sponge';

    const hasGenie3Data = meta?.[genie3Key] && Object.keys(meta[genie3Key]).length > 0;
    const hasSpongeData = meta?.[spongeKey] && Object.keys(meta[spongeKey]).length > 0;

    if (!hasGenie3Data && !hasSpongeData) {
      console.warn(
        `[RegulatoryScores] No ${displayMode} data available for ${compare ? 'compare' : 'main'} view. Available keys:`,
        Object.keys(meta || {}).filter(k => k.includes('regulatory'))
      );
    }

    console.log(
      `[RegulatoryScores] Display mode changed to '${displayMode}' for ${compare ? 'compare' : 'main'} view`
    );
  }

  public getRegulatoryActionColumns(networkType: 'genie3' | 'sponge', compare: boolean): string[] {
    const meta = compare ? this.metaCompare : this.meta;
    const scoreNames = meta?.['grn_score_names'];
    if (!Array.isArray(scoreNames)) {
      return [];
    }
    return scoreNames.filter((name: string) => name.endsWith(`_${networkType}`));
  }

  public getRegulatoryTableData(networkType: 'genie3' | 'sponge', compare: boolean): TableData {
    const displayMode = compare ? this.regulatoryScoreDisplayModeCompare : this.regulatoryScoreDisplayMode;
    const meta = compare ? this.metaCompare : this.meta;
    const key =
      displayMode === 'raw'
        ? `global_regulatory_scores_${networkType}`
        : displayMode === 'moranI'
          ? `global_regulatory_moranI_${networkType}`
          : `global_regulatory_gearyC_${networkType}`;

    const table = meta?.[key];
    if (!table || typeof table !== 'object') {
      return {};
    }
    return table as TableData;
  }

  public getRegulatoryEmptyMessage(networkType: 'genie3' | 'sponge', compare: boolean): string {
    const displayMode = compare ? this.regulatoryScoreDisplayModeCompare : this.regulatoryScoreDisplayMode;
    if (displayMode === 'moranI') {
      return `No ${networkType === 'genie3' ? 'Genie3' : 'Sponge'} Moran's I regulatory scores loaded`;
    }
    if (displayMode === 'gearyC') {
      return `No ${networkType === 'genie3' ? 'Genie3' : 'Sponge'} Geary's C regulatory scores loaded`;
    }
    return `No ${networkType === 'genie3' ? 'Genie3' : 'Sponge'} global regulatory scores loaded`;
  }
  // isArray/isNestedArray/isNumberLike/toNumberLike/isNumericArray/isPrimitive/getArrayStats/
  // formatValue/toggleExpand/isExpanded/trackByIndex/asArrayRow/looksLikeArrayString/
  // toJsonCompact/shouldShowProperty/groupedCellProperties/computeGroupedProperties/
  // isPlainObject/objectKeyCount/objectEntries/prettyKey/dictId/expandedProps/hiddenPropKeys all
  // moved into CellInfoPanelComponent along with the ng-template they backed — see
  // cell-info-panel/cell-info-panel.component.ts. `label()` below is kept here too (small,
  // duplicated, pure) since this component's own renderLegend still needs it.

  label(key: string): string {
    const translated = this.translationService.translateSync(key);
    if (translated !== key) return translated;
    return key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
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

  public selectDatasetTutorial(): void {
    const tour = new Shepherd.Tour({
      useModalOverlay: true,
      defaultStepOptions: {
        classes: 'shepherd-theme-custom'
      }
    });

    tour.addStep({
      id: 'dataset-dropdown',
      attachTo: { element: '#dataset-dropdown-main', on: 'left' },
      text: 'This dropdown allows you to select the dataset you want to explore.  It shows:' +
      ' <ul><li>All of your uploaded and the builtin datasets</li>' +
      ' <li>Tangram datasets will allow you to also select to visualize the mapped data as well (under Tangram)</li></ul>',
      buttons: [{ text: "Done", action: tour.complete }]
    });

    tour.start();
  }

  public globalRegulatoryScoresTutorial(): void {
    const tour = new Shepherd.Tour({
      useModalOverlay: true,
      defaultStepOptions: {
        classes: 'shepherd-theme-custom'
      }
    });

    tour.addStep({
      id: 'global-regulatory-scores',
      attachTo: { element: '#global-regulatory-scores-info', on: 'left' },
      text: 'The global regulatory scores summarize overall regulatory activity per gene set by aggregating the scores over all cells. They can be used to identify genesets with overall high regulatory activity or to find general candidate genesets to prioritize for further analysis.',
      buttons: [{ text: "Done", action: tour.complete }]
    });

    tour.start();
  }

  private renderLegend(containerName: string): void {
    const viewVariablesToUpdate = this.getViewVariablesToUpdate(containerName);
    const svg = viewVariablesToUpdate.svg;
    // svg now comes from the mainView/compareView @ViewChild ref (see getViewVariablesToUpdate),
    // so it's typed as possibly undefined if the child hasn't been created yet — matches the
    // original code's implicit assumption that it's always set by the time this runs.
    if (!svg) return;

    svg.selectAll(`.${viewVariablesToUpdate.legendContainerName}`).remove();

    if (viewVariablesToUpdate.getLegendType() === 'continuous') {
      const [min, max] = viewVariablesToUpdate.getLegendDomain() as number[] || [0, 1];
      const legendX = 20;
      const legendY = 20;
      const width = 240;
      const height = 24;
      const fontSize = 18;
      const padding = 14;

      // Use standard <defs>
      const defs = svg.select('defs').empty()
        ? svg.append('defs')
        : svg.select('defs');

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

      const legendG = svg
        .append('g')
        .attr('class', viewVariablesToUpdate.legendContainerName)
        .attr('transform', `translate(${legendX},${legendY})`);

      const titleText = this.translationService.translateSync(viewVariablesToUpdate.view);
      // fallback if translation returns empty
      const legendTitle = titleText && String(titleText).trim() ? titleText : this.label(viewVariablesToUpdate.view);

      // measure sizes
      const tempSvg = svg.append('g').style('opacity', 0);
      const titleWidth =
        tempSvg
          .append('text')
          .text(legendTitle)
          .style('font-size', `${fontSize}px`)
          .style('font-weight', 'bold')
          .node()
          ?.getBBox().width || 0;

      const minText = (min ?? 0).toFixed(2);
      const minWidth =
        tempSvg
          .append('text')
          .text(minText)
          .style('font-size', `${fontSize - 2}px`)
          .node()
          ?.getBBox().width || 0;

      const maxText = (max ?? 0).toFixed(2);
      const maxWidth =
        tempSvg
          .append('text')
          .text(maxText)
          .style('font-size', `${fontSize - 2}px`)
          .node()
          ?.getBBox().width || 0;

      tempSvg.remove();

      const textHeight = fontSize * 1.2;
      const contentWidth = Math.max(width, titleWidth, minWidth + maxWidth + 20);
      const bgWidth = contentWidth + padding * 2;
      const bgHeight = padding + textHeight + 6 + height + 4 + textHeight + padding;

      legendG
        .append('rect')
        .attr('x', 0)
        .attr('y', 0)
        .attr('width', bgWidth)
        .attr('height', bgHeight)
        .style('fill', 'rgba(255, 255, 255, 0.95)')
        .attr('stroke', '#d0d0d0')
        .attr('stroke-width', 1.2)
        .attr('rx', 8)
        .style('filter', 'drop-shadow(0 2px 6px rgba(0,0,0,0.08))');

      const titleY = padding + fontSize;
      legendG
        .append('text')
        .attr('x', bgWidth / 2)
        .attr('y', titleY)
        .attr('text-anchor', 'middle')
        .style('font-size', `${fontSize}px`)
        .style('font-weight', 'bold')
        .style('fill', '#222')
        .text(legendTitle);

      const barX = (bgWidth - contentWidth) / 2;
      const barY = titleY + 8;
      legendG
        .append('rect')
        .attr('x', barX)
        .attr('y', barY)
        .attr('width', contentWidth)
        .attr('height', height)
        .style('fill', `url(#${viewVariablesToUpdate.legendGradientName})`)
        .attr('stroke', '#ccc')
        .attr('stroke-width', 1)
        .attr('rx', 4);

      const labelsY = barY + height + fontSize + 2;
      // Min label
      legendG
        .append('text')
        .attr('x', barX)
        .attr('y', labelsY)
        .attr('text-anchor', 'start')
        .style('font-size', `${fontSize - 2}px`)
        .style('fill', '#444')
        .text(minText);

      // Max label
      legendG
        .append('text')
        .attr('x', barX + contentWidth)
        .attr('y', labelsY)
        .attr('text-anchor', 'end')
        .style('font-size', `${fontSize - 2}px`)
        .style('fill', '#444')
        .text(maxText);

    } else {
      // Categorical legend
      const categories = viewVariablesToUpdate.getLegendDomain() as string[] || [];
      categories.sort();
      const legendX = 20;
      const legendY = 20;
      const itemHeight = 32;
      const rectHeight = 16;
      const rectWidth = 24;
      const fontSize = 18;
      const padding = 14;

      // Measure using svg
      const tempSvg = svg.append('g').style('opacity', 0);
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

      const itemWidth = Math.max(160, maxTextWidth + rectWidth + 14, titleWidth);
      const backgroundWidth = itemWidth + padding * 2;
      const titleHeight = fontSize * 1.3;
      const categoryItemsHeight = categories.length * itemHeight;
      const backgroundHeight = padding + titleHeight + 8 + categoryItemsHeight + padding;

      const legendG = svg
        .append('g')
        .attr('class', viewVariablesToUpdate.legendContainerName)
        .attr('transform', `translate(${legendX},${legendY})`);

      // Background rectangle for categorical legend
      legendG
        .append('rect')
        .attr('x', 0)
        .attr('y', 0)
        .attr('width', backgroundWidth)
        .attr('height', backgroundHeight)
        .style('fill', 'rgba(255, 255, 255, 0.95)')
        .attr('stroke', '#d0d0d0')
        .attr('stroke-width', 1.2)
        .attr('rx', 8)
        .style('filter', 'drop-shadow(0 2px 6px rgba(0,0,0,0.08))');

      // Title (categorical)
      const titleY = padding + fontSize;
      legendG
        .append('text')
        .attr('x', backgroundWidth / 2)
        .attr('y', titleY)
        .attr('text-anchor', 'middle')
        .style('font-size', `${fontSize}px`)
        .style('font-weight', 'bold')
        .style('fill', '#222')
        .text(legendTitleCat);

      const itemsStartY = padding + titleHeight + 8;
      categories.forEach((cat, i) => {
        const yPosition = itemsStartY + i * itemHeight;
        const legendItem = legendG.append('g').attr('transform', `translate(${padding}, ${yPosition})`);
        const rectY = (itemHeight - rectHeight) / 2;
        legendItem
          .append('rect')
          .attr('x', 0)
          .attr('y', rectY)
          .attr('width', rectWidth)
          .attr('height', rectHeight)
          .style('fill', viewVariablesToUpdate.ordinal(cat))
          .attr('stroke', '#333')
          .attr('stroke-width', 0.5)
          .attr('rx', 3);
        const textY = rectY + rectHeight / 2;
        legendItem
          .append('text')
          .attr('x', rectWidth + 10)
          .attr('y', textY)
          .attr('dy', '0.35em')
          .style('font-size', `${fontSize - 2}px`)
          .style('fill', '#333')
          .text(cat);
        legendItem.append('title').text(cat);
      });
    }
  }

  /**
   * Load TF graph data from the backend API
   */
  public loadTfGraph(compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;

    if (!dataset?.id) {
      console.warn('No dataset selected');
      return;
    }

    const loading = compare ? 'tfGraphLoadingCompare' : 'tfGraphLoading';
    this[loading] = true;

    // Build the API URL
    const apiUrl = `/api/tf_graph/${dataset.id}`;
    const fullUrl = `${this.sessionService.apiUrl}${apiUrl}`;

    console.log('[TF Graph] Loading from:', fullUrl);
    console.log('[TF Graph] Loading flag set for', compare ? 'compare' : 'main');

    // Ensure container exists in DOM
    const containerId = compare ? '#tf_gene_net_container_compare' : '#tf_gene_net_container';
    const container = document.querySelector(containerId);
    console.log('[TF Graph] Container exists:', !!container, 'ID:', containerId);

    // Load the JSON file using fetch like GeoJSON does
    fetch(fullUrl, { credentials: 'include' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`HTTP Error: ${response.status} ${response.statusText}`);
        }
        return response.json();
      })
      .then((data: { nodes: any[], links: any[] }) => {
        console.log('✓ Loaded TF graph data:', data);
        this.visualizeTfGraph(data, compare);
        this[loading] = false;
      })
      .catch((err) => {
        console.error('✗ Failed to load TF graph:', err);
        this[loading] = false;
      });
  }

  /**
   * Visualize TF graph using D3
   */
  private visualizeTfGraph(graphData: { nodes: any[], links: any[] }, compare: boolean = false): void {
    const containerId = compare ? '#tf_gene_net_container_compare' : '#tf_gene_net_container';

    console.log('[visualizeTfGraph] Selecting container:', containerId);
    const container = d3.select(containerId);
    console.log('[visualizeTfGraph] Container selected:', container.empty() ? 'EMPTY' : 'EXISTS');

    // Clear previous content
    container.html('');

    if (!graphData.nodes || !graphData.links || graphData.nodes.length === 0) {
      console.warn('[visualizeTfGraph] No graph data available');
      container.append('p').text('No TF graph data available');
      return;
    }

    console.log('[visualizeTfGraph] Nodes:', graphData.nodes.length, 'Links:', graphData.links.length);

    // SVG dimensions
    const width = 1000;
    const height = 600;

    // Create SVG
    const svg = d3
      .select(containerId)
      .append('svg')
      .attr('width', '100%')
      .attr('height', height)
      .style('background-color', '#f8f9fa')
      .style('border', '1px solid #ddd');

    console.log('[visualizeTfGraph] SVG created and appended');

    // Add zoom behavior
    const g = svg.append('g');
    const zoom = d3.zoom<SVGSVGElement, any>().on('zoom', (event) => {
      g.attr('transform', event.transform);
    });
    svg.call(zoom);

    // Create map from node id to node object for efficient lookup
    const nodeMap = new Map(graphData.nodes.map((node: any) => [node.id, node]));

    // Transform links to use node objects instead of just ids
    const links = graphData.links.map((link: any) => ({
      source: nodeMap.get(link.source),
      target: nodeMap.get(link.target),
      edge_color: link.edge_color || '#999',
      edge_width: link.edge_width || 1,
      edge_dash: link.edge_dash || '0'
    })).filter((link: any) => link.source && link.target);

    console.log('[visualizeTfGraph] Filtered links:', links.length);

    // Force simulation
    const simulation = d3
      .forceSimulation(graphData.nodes)
      .force(
        'link',
        d3
          .forceLink(links)
          .id((d: any) => d.id)
          .distance(50)
          .strength(0.4)
      )
      .force('charge', d3.forceManyBody().strength(-500))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(30));

    // Draw links
    const link = g
      .append('g')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.3)
      .selectAll('line')
      .data(links)
      .enter()
      .append('line')
      .attr('stroke', (d: any) => d.edge_color)
      .attr('stroke-width', (d: any) => d.edge_width)
      .style('stroke-dasharray', (d: any) => d.edge_dash);

    // Draw nodes
    const node = g
      .append('g')
      .attr('stroke', '#fff')
      .attr('stroke-width', 2)
      .selectAll('circle')
      .data(graphData.nodes)
      .enter()
      .append('circle')
      .attr('r', (d: any) => d.type === 'from' ? 8 : 6)
      .attr('fill', (d: any) => d.node_fill || '#999')
      .attr('opacity', 0.9)
      .call(
        d3
          .drag<SVGCircleElement, any>()
          .on('start', (event: any, d: any) => {
            if (!event.active) simulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
          })
          .on('drag', (event: any, d: any) => {
            d.fx = event.x;
            d.fy = event.y;
          })
          .on('end', (event: any, d: any) => {
            if (!event.active) simulation.alphaTarget(0);
            d.fx = null;
            d.fy = null;
          })
      );

    // Add node labels
    const labels = g
      .append('g')
      .selectAll('text')
      .data(graphData.nodes)
      .enter()
      .append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '.35em')
      .attr('font-size', (d: any) => d.type === 'from' ? '11px' : '10px')
      .attr('font-weight', (d: any) => d.priorTF ? 'bold' : 'normal')
      .attr('fill', '#333')
      .attr('pointer-events', 'none')
      .text((d: any) => d.name.length > 10 ? d.name.substring(0, 10) + '...' : d.name);

    // Add tooltips on hover
    node.append('title').text((d: any) => `${d.name}${d.priorTF ? ' (Prior TF)' : ''}`);

    // Update positions on simulation tick
    simulation.on('tick', () => {
      link
        .attr('x1', (d: any) => d.source.x)
        .attr('y1', (d: any) => d.source.y)
        .attr('x2', (d: any) => d.target.x)
        .attr('y2', (d: any) => d.target.y);

      node
        .attr('cx', (d: any) => d.x)
        .attr('cy', (d: any) => d.y);

      labels
        .attr('x', (d: any) => d.x)
        .attr('y', (d: any) => d.y);
    });

    // Add legend
    this.addTfGraphLegend(svg, width, height);
  }

  /**
   * Add legend to TF graph visualization
   */
  private addTfGraphLegend(svg: any, width: number, height: number): void {
    const legendX = 20;
    const legendY = 20;

    const legend = svg.append('g')
      .attr('transform', `translate(${legendX}, ${legendY})`);

    // Background
    legend.append('rect')
      .attr('width', 180)
      .attr('height', 100)
      .attr('fill', 'white')
      .attr('stroke', '#ccc')
      .attr('rx', 4);

    // Title
    legend.append('text')
      .attr('x', 10)
      .attr('y', 20)
      .attr('font-weight', 'bold')
      .attr('font-size', '12px')
      .text('Legend');

    // TF nodes
    legend.append('circle')
      .attr('cx', 20)
      .attr('cy', 40)
      .attr('r', 5)
      .attr('fill', '#1f77b4');

    legend.append('text')
      .attr('x', 35)
      .attr('y', 45)
      .attr('font-size', '11px')
      .text('Regulator');

    // Prior TF indicator
    legend.append('circle')
      .attr('cx', 20)
      .attr('cy', 60)
      .attr('r', 5)
      .attr('fill', '#ff7f0e');

    legend.append('text')
      .attr('x', 35)
      .attr('y', 65)
      .attr('font-size', '11px')
      .text('Prior TF');

    // Target nodes
    legend.append('circle')
      .attr('cx', 20)
      .attr('cy', 80)
      .attr('r', 4)
      .attr('fill', '#cccccc');

    legend.append('text')
      .attr('x', 35)
      .attr('y', 85)
      .attr('font-size', '11px')
      .text('Target');
  }

  /**
   * Load pre-computed GRN graph (from dataset metadata)
   * Loads all 3 filter modes in parallel and stores them for display
   */
  public loadPrecomputedGrnGraph(compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;

    if (!dataset?.id) {
      console.warn('[loadPrecomputedGrnGraph] No dataset selected');
      return;
    }

    // Get the GRN evaluation name from dataset metadata or use default
    const evalName = dataset.grn_evaluation_name || 'GRN_Evaluation';
    console.log('[loadPrecomputedGrnGraph] Loading for dataset:', dataset.id, 'evaluation:', evalName);

    // Load all 3 filter modes in parallel
    const filterModes: Array<'prior' | 'extended' | 'full'> = ['prior', 'extended', 'full'];
    const self = this as any;

    filterModes.forEach((filterMode) => {
      const loadingKey = compare
        ? `isLoadingGrnGraph${filterMode.charAt(0).toUpperCase() + filterMode.slice(1)}Compare`
        : `isLoadingGrnGraph${filterMode.charAt(0).toUpperCase() + filterMode.slice(1)}`;
      const dataKey = compare
        ? `grnGraphData${filterMode.charAt(0).toUpperCase() + filterMode.slice(1)}Compare`
        : `grnGraphData${filterMode.charAt(0).toUpperCase() + filterMode.slice(1)}`;

      self[loadingKey] = true;

      const url = `/api/grn_evaluation/${dataset.id}/${evalName}/graph?filter_mode=${filterMode}`;
      console.log(`[loadPrecomputedGrnGraph] Fetching ${filterMode}: ${url}`);

      fetch(url, { credentials: 'include' })
        .then((response) => {
          if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${response.statusText}`);
          }
          return response.json();
        })
        .then((data) => {
          console.log(`[loadPrecomputedGrnGraph] Loaded ${filterMode} data:`, data?.nodes?.length, 'nodes,', data?.links?.length, 'links');
          self[dataKey] = data;
          self[loadingKey] = false;

          // Check if this is the currently selected filter and auto-visualize only then
          const currentFilterKey = compare ? 'grnGraphFilterCompare' : 'grnGraphFilter';
          if (self[currentFilterKey] === filterMode) {
            console.log(`[loadPrecomputedGrnGraph] Auto-visualizing ${filterMode} (current selection)`);
            // Use setTimeout to ensure DOM is ready
            setTimeout(() => {
              this.visualizeGrnGraph(data, compare, filterMode);
            }, 50);
          } else {
            console.log(`[loadPrecomputedGrnGraph] Cached ${filterMode} data, not visualizing (current filter is ${self[currentFilterKey]})`);
          }
        })
        .catch((err) => {
          console.error(`[loadPrecomputedGrnGraph] Error loading ${filterMode}:`, err);
          self[loadingKey] = false;
        });
    });
  }

  /**
   * Load pre-computed GRN plots (peak_plots and motif_plots directories)
   */
  public loadPrecomputedGrnPlots(compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;

    if (!dataset?.id) {
      console.warn('[loadPrecomputedGrnPlots] No dataset selected');
      return;
    }

    // Get the GRN evaluation name from dataset metadata or use default
    const evalName = dataset.grn_evaluation_name || 'GRN_Evaluation';
    console.log('[loadPrecomputedGrnPlots] Loading for dataset:', dataset.id, 'evaluation:', evalName);

    const plotsKey = compare ? 'grnPlotsCompare' : 'grnPlots';
    const isLoadingKey = compare ? 'isLoadingGrnPlotsCompare' : 'isLoadingGrnPlots';

    const self = this as any;
    self[isLoadingKey] = true;

    const url = `/api/grn_evaluation/${dataset.id}/${evalName}/plots`;
    console.log(`[loadPrecomputedGrnPlots] Fetching from: ${url}`);

    fetch(url, { credentials: 'include' })
      .then((response) => {
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}: ${response.statusText}`);
        }
        return response.json();
      })
      .then((data) => {
        console.log('[loadPrecomputedGrnPlots] Loaded plots:', data);
        self[plotsKey] = data;
        self[isLoadingKey] = false;
      })
      .catch((err) => {
        console.error('[loadPrecomputedGrnPlots] Error loading plots:', err);
        self[isLoadingKey] = false;
        self[plotsKey] = {}; // Empty object on error
      });
  }

  computeOnDemandGrnEvaluation(compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) {
      compare ? this.onDemandGrnErrorCompare = 'No dataset selected' : this.onDemandGrnError = 'No dataset selected';
      return;
    }

    const params = compare ? this.onDemandGrnParamsCompare : this.onDemandGrnParams;
    compare ? this.isComputingOnDemandGrnCompare = true : this.isComputingOnDemandGrn = true;
    compare ? this.onDemandGrnErrorCompare = '' : this.onDemandGrnError = '';

    const payload = {
      dataset_id: dataset.id,
      grn_evaluation_obs_key: params.obsKey,
      grn_evaluation_cluster: params.cluster,
      grn_evaluation_name: params.name,
      grn_evaluation_gene_set: params.geneSet,
      grn_evaluation_use_prior_grn: params.usePriorGrn,
    };

    this.http.post('/api/on_demand_grn_evaluation', payload).subscribe({
      next: (response: any) => {
        compare ? this.isComputingOnDemandGrnCompare = false : this.isComputingOnDemandGrn = false;
        alert('GRN Evaluation computed successfully');
        if (compare) {
          this.onDemandOutputPathCompare = response.output_path;
        } else {
          this.onDemandOutputPath = response.output_path;
        }

      },
      error: (err) => {
        const message = err?.error?.detail || err.message || 'Error computing GRN Evaluation';
        if (compare) {
          this.isComputingOnDemandGrnCompare = false;
          this.onDemandGrnErrorCompare = message;
        } else {
          this.isComputingOnDemandGrn = false;
          this.onDemandGrnError = message;
        }
      }
    });
  }

  public loadOnDemandGrnRuns(compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) return;

    const isLoadingKey = compare ? 'isLoadingGrnRunsCompare' : 'isLoadingGrnRuns';
    const runsKey = compare ? 'onDemandGrnRunsCompare' : 'onDemandGrnRuns';
    const errorKey = compare ? 'grnTabErrorCompare' : 'grnTabError';

    this[isLoadingKey] = true;
    this[errorKey] = null;

    this.http.get<{ runs: string[] }>(`/api/grn_evaluation/${dataset.id}/available_runs`)
      .subscribe({
        next: (response) => {
          this[runsKey] = response.runs;
          this[isLoadingKey] = false;
        },
        error: (err) => {
          this[errorKey] = 'Failed to load available GRN runs';
          this[isLoadingKey] = false;
        }
      });
  }

  /**
   * Handle when user selects a GRN run
   */
  public onSelectGrnRun(runName: string, compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) return;

    const selectedRunKey = compare ? 'selectedOnDemandRunCompare' : 'selectedOnDemandRun';
    const errorKey = compare ? 'grnTabErrorCompare' : 'grnTabError';

    this[selectedRunKey] = runName;
    this[errorKey] = null;

    // Load graph, plots, and stats in parallel
    this.loadGrnGraph(runName, compare);
    this.loadGrnStats(runName, compare);
    this.loadGrnPlots(runName, compare);
  }

  /**
   * Load all 3 GRN graphs (prior, extended, full) and display selected one
   */
  private loadGrnGraph(runName: string, compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) {
      console.warn('[loadGrnGraph] No dataset selected');
      return;
    }

    console.log('[loadGrnGraph] Loading GRN graphs for run:', runName, 'Dataset:', dataset.id, 'Compare:', compare);

    // Load PRIOR
    this.isLoadingGrnGraphPrior = !compare ? true : this.isLoadingGrnGraphPrior;
    if (compare) this.isLoadingGrnGraphPriorCompare = true;
    this.http.get<any>(
      `/api/grn_evaluation/${dataset.id}/${runName}/graph`,
      { params: { filter_mode: 'prior' } }
    ).subscribe({
      next: (graphData) => {
        console.log('[loadGrnGraph] Loaded PRIOR data. Nodes:', graphData?.nodes?.length, 'Links:', graphData?.links?.length);
        if (compare) {
          this.grnGraphDataPriorCompare = graphData;
          this.isLoadingGrnGraphPriorCompare = false;
        } else {
          this.grnGraphDataPrior = graphData;
          this.isLoadingGrnGraphPrior = false;
        }
        this.visualizeGrnGraphForFilter('prior', compare);
      },
      error: (err) => {
        console.error('[loadGrnGraph] Error loading PRIOR:', err);
        if (compare) this.isLoadingGrnGraphPriorCompare = false;
        else this.isLoadingGrnGraphPrior = false;
      }
    });

    // Load EXTENDED
    this.isLoadingGrnGraphExtended = !compare ? true : this.isLoadingGrnGraphExtended;
    if (compare) this.isLoadingGrnGraphExtendedCompare = true;
    this.http.get<any>(
      `/api/grn_evaluation/${dataset.id}/${runName}/graph`,
      { params: { filter_mode: 'extended' } }
    ).subscribe({
      next: (graphData) => {
        console.log('[loadGrnGraph] Loaded EXTENDED data. Nodes:', graphData?.nodes?.length, 'Links:', graphData?.links?.length);
        if (compare) {
          this.grnGraphDataExtendedCompare = graphData;
          this.isLoadingGrnGraphExtendedCompare = false;
        } else {
          this.grnGraphDataExtended = graphData;
          this.isLoadingGrnGraphExtended = false;
        }
        this.visualizeGrnGraphForFilter('extended', compare);
      },
      error: (err) => {
        console.error('[loadGrnGraph] Error loading EXTENDED:', err);
        if (compare) this.isLoadingGrnGraphExtendedCompare = false;
        else this.isLoadingGrnGraphExtended = false;
      }
    });

    // Load FULL
    this.isLoadingGrnGraphFull = !compare ? true : this.isLoadingGrnGraphFull;
    if (compare) this.isLoadingGrnGraphFullCompare = true;
    this.http.get<any>(
      `/api/grn_evaluation/${dataset.id}/${runName}/graph`,
      { params: { filter_mode: 'full' } }
    ).subscribe({
      next: (graphData) => {
        console.log('[loadGrnGraph] Loaded FULL data. Nodes:', graphData?.nodes?.length, 'Links:', graphData?.links?.length);
        if (compare) {
          this.grnGraphDataFullCompare = graphData;
          this.isLoadingGrnGraphFullCompare = false;
        } else {
          this.grnGraphDataFull = graphData;
          this.isLoadingGrnGraphFull = false;
        }
        this.visualizeGrnGraphForFilter('full', compare);
      },
      error: (err) => {
        console.error('[loadGrnGraph] Error loading FULL:', err);
        if (compare) this.isLoadingGrnGraphFullCompare = false;
        else this.isLoadingGrnGraphFull = false;
      }
    });
  }

  /**
   * Visualize the currently selected GRN graph filter
   */
  private visualizeGrnGraphForFilter(filterMode: 'prior' | 'extended' | 'full', compare: boolean = false): void {
    const currentFilter = compare ? this.grnGraphFilterCompare : this.grnGraphFilter;

    console.log('[visualizeGrnGraphForFilter] Mode:', filterMode, 'CurrentFilter:', currentFilter, 'Match:', filterMode === currentFilter);

    // Only visualize if this is the currently selected filter
    if (filterMode === currentFilter) {
      let graphData: any = null;

      if (compare) {
        if (filterMode === 'prior') graphData = this.grnGraphDataPriorCompare;
        else if (filterMode === 'extended') graphData = this.grnGraphDataExtendedCompare;
        else if (filterMode === 'full') graphData = this.grnGraphDataFullCompare;
      } else {
        if (filterMode === 'prior') graphData = this.grnGraphDataPrior;
        else if (filterMode === 'extended') graphData = this.grnGraphDataExtended;
        else if (filterMode === 'full') graphData = this.grnGraphDataFull;
      }

      console.log('[visualizeGrnGraphForFilter] Graph data exists:', !!graphData, 'Nodes:', graphData?.nodes?.length);

      if (graphData) {
        this.visualizeGrnGraph(graphData, compare);
      } else {
        console.warn('[visualizeGrnGraphForFilter] No graph data for mode:', filterMode);
      }
    }
  }

  /**
   * Change graph filter mode and switch visualization
   */
  public onChangeGrnGraphFilter(newFilter: 'prior' | 'extended' | 'full', compare: boolean = false): void {
    console.log('[onChangeGrnGraphFilter] Changing filter to:', newFilter, 'Compare:', compare);

    if (compare) {
      this.grnGraphFilterCompare = newFilter;
    } else {
      this.grnGraphFilter = newFilter;
    }

    // Get the data for the selected filter
    let graphData: any = null;

    if (compare) {
      if (newFilter === 'prior') graphData = this.grnGraphDataPriorCompare;
      else if (newFilter === 'extended') graphData = this.grnGraphDataExtendedCompare;
      else if (newFilter === 'full') graphData = this.grnGraphDataFullCompare;
    } else {
      if (newFilter === 'prior') graphData = this.grnGraphDataPrior;
      else if (newFilter === 'extended') graphData = this.grnGraphDataExtended;
      else if (newFilter === 'full') graphData = this.grnGraphDataFull;
    }

    console.log('[onChangeGrnGraphFilter] Data exists for', newFilter + ':', !!graphData);

    if (graphData) {
      console.log('[onChangeGrnGraphFilter] Rendering immediately for filter:', newFilter);
      this.visualizeGrnGraph(graphData, compare, newFilter);
    } else {
      // Otherwise load from backend
      const runName = compare ? this.selectedOnDemandRunCompare : this.selectedOnDemandRun;
    }
  }

  /**
   * Load available plot images
   */
  private loadGrnPlots(runName: string, compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) return;

    const isLoadingKey = compare ? 'isLoadingGrnPlotsCompare' : 'isLoadingGrnPlots';
    const plotsKey = compare ? 'grnPlotsCompare' : 'grnPlots';

    this[isLoadingKey] = true;

    this.http.get<{ [dir: string]: string[] }>(
      `/api/grn_evaluation/${dataset.id}/${runName}/plots`
    ).subscribe({
      next: (plots) => {
        this[plotsKey] = plots;
        this[isLoadingKey] = false;
      },
      error: (err) => {
        this[plotsKey] = {};
        this[isLoadingKey] = false;
      }
    });
  }

  /**
   * Load peak and motif stats
   */
  private loadGrnStats(runName: string, compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) return;

    const isLoadingKey = compare ? 'isLoadingGrnStatsCompare' : 'isLoadingGrnStats';
    const peakStatsKey = compare ? 'grnPeakStatsCompare' : 'grnPeakStats';
    const motifStatsKey = compare ? 'grnMotifStatsCompare' : 'grnMotifStats';

    this[isLoadingKey] = true;

    // Load peak stats
    this.http.get<any>(
      `/api/grn_evaluation/${dataset.id}/${runName}/stats/peak`
    ).subscribe({
      next: (stats) => {
        this[peakStatsKey] = stats;
      },
      error: () => {
        this[peakStatsKey] = null;
      }
    });

    // Load motif stats
    this.http.get<any>(
      `/api/grn_evaluation/${dataset.id}/${runName}/stats/motif`
    ).subscribe({
      next: (stats) => {
        this[motifStatsKey] = stats;
        this[isLoadingKey] = false;
      },
      error: () => {
        this[motifStatsKey] = null;
        this[isLoadingKey] = false;
      }
    });
  }

  /**
   * Visualize GRN graph using D3 (reuse TF Gene-Net logic)
   */
  private visualizeGrnGraph(graphData: any, compare: boolean = false, filterMode?: 'prior' | 'extended' | 'full'): void {
    // If filterMode not provided, get it from component state
    if (!filterMode) {
      const filterKey = compare ? 'grnGraphFilterCompare' : 'grnGraphFilter';
      filterMode = this[filterKey];
    }

    // Determine container ID based on filter mode and compare status
    let containerId = '#tf_gene_net_container_' + filterMode;
    if (compare) {
      containerId = '#tf_gene_net_container_' + filterMode + '_compare';
    }

    console.log('[visualizeGrnGraph] Starting visualization');
    console.log('[visualizeGrnGraph] Filter mode:', filterMode);
    console.log('[visualizeGrnGraph] Container ID:', containerId);
    console.log('[visualizeGrnGraph] Compare:', compare);
    console.log('[visualizeGrnGraph] Graph data:', graphData);

    // Check if container exists in DOM
    const containerElement = document.querySelector(containerId);
    console.log('[visualizeGrnGraph] Container element found in DOM:', !!containerElement);
    if (!containerElement) {
      console.error('[visualizeGrnGraph] Container not found in DOM! ID:', containerId);
      return;
    }

    // Clear previous content
    d3.select(containerId).html('');
    console.log('[visualizeGrnGraph] Cleared previous content');

    if (!graphData || !graphData.nodes || !graphData.links || graphData.nodes.length === 0) {
      console.warn('[visualizeGrnGraph] No graph data available for mode:', filterMode);
      d3.select(containerId).append('p').text('No GRN graph data available');
      return;
    }

    console.log('[visualizeGrnGraph] Rendering graph with', graphData.nodes.length, 'nodes and', graphData.links.length, 'links');

    // SVG dimensions
    const width = 1000;
    const height = 600;

    // Create SVG
    const svg = d3
      .select(containerId)
      .append('svg')
      .attr('width', '100%')
      .attr('height', height)
      .style('background-color', '#f8f9fa')
      .style('border', '1px solid #ddd');

    // Add zoom behavior
    const g = svg.append('g');
    const zoom = d3.zoom<SVGSVGElement, any>().on('zoom', (event) => {
      g.attr('transform', event.transform);
    });
    svg.call(zoom);

    // Create map from node id to node object
    const nodeMap = new Map(graphData.nodes.map((node: any) => [node.id, node]));

    // Transform links to use node objects
    const links = graphData.links.map((link: any) => ({
      source: nodeMap.get(link.source),
      target: nodeMap.get(link.target),
      edge_color: link.edge_color || '#999',
      edge_width: link.edge_width || 1,
      edge_dash: link.edge_dash || '0'
    })).filter((link: any) => link.source && link.target);

    // Force simulation
    const simulation = d3
      .forceSimulation(graphData.nodes)
      .force(
        'link',
        d3
          .forceLink(links)
          .id((d: any) => d.id)
          .distance(50)
          .strength(0.4)
      )
      .force('charge', d3.forceManyBody().strength(-500))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(30));

    // Draw links
    const link = g
      .append('g')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.3)
      .selectAll('line')
      .data(links)
      .enter()
      .append('line')
      .attr('stroke', (d: any) => d.edge_color)
      .attr('stroke-width', (d: any) => d.edge_width)
      .style('stroke-dasharray', (d: any) => d.edge_dash);

    // Draw nodes
    const node = g
      .append('g')
      .attr('stroke', '#fff')
      .attr('stroke-width', 2)
      .selectAll('circle')
      .data(graphData.nodes)
      .enter()
      .append('circle')
      .attr('r', (d: any) => d.type === 'from' ? 8 : 6)
      .attr('fill', (d: any) => d.node_fill || '#999')
      .attr('opacity', 0.9)
      .call(
        d3
          .drag<SVGCircleElement, any>()
          .on('start', (event: any, d: any) => {
            if (!event.active) simulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
          })
          .on('drag', (event: any, d: any) => {
            d.fx = event.x;
            d.fy = event.y;
          })
          .on('end', (event: any, d: any) => {
            if (!event.active) simulation.alphaTarget(0);
            d.fx = null;
            d.fy = null;
          })
      );

    // Add node labels
    const labels = g
      .append('g')
      .selectAll('text')
      .data(graphData.nodes)
      .enter()
      .append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '.35em')
      .attr('font-size', (d: any) => d.type === 'from' ? '11px' : '10px')
      .attr('font-weight', (d: any) => d.priorTF ? 'bold' : 'normal')
      .attr('fill', '#333')
      .attr('pointer-events', 'none')
      .text((d: any) => d.name.length > 10 ? d.name.substring(0, 10) + '...' : d.name);

    // Add tooltips
    node.append('title').text((d: any) => `${d.name}${d.priorTF ? ' (Prior TF)' : ''}`);

    // Update positions on tick
    simulation.on('tick', () => {
      link
        .attr('x1', (d: any) => d.source.x)
        .attr('y1', (d: any) => d.source.y)
        .attr('x2', (d: any) => d.target.x)
        .attr('y2', (d: any) => d.target.y);

      node
        .attr('cx', (d: any) => d.x)
        .attr('cy', (d: any) => d.y);

      labels
        .attr('x', (d: any) => d.x)
        .attr('y', (d: any) => d.y);
    });
  }

  /**
   * Get image URL for a plot
   */
  public getGrnImageUrl(imageName: string, compare: boolean = false): string {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    if (!dataset) return '';

    // For on-demand runs, use the selected run name
    const selectedRun = compare ? this.selectedOnDemandRunCompare : this.selectedOnDemandRun;
    if (selectedRun) {
      return `/api/grn_evaluation/${dataset.id}/${selectedRun}/image/${imageName}`;
    }

    // For pre-computed plots, use the dataset's GRN evaluation name
    const evalName = dataset.grn_evaluation_name || 'GRN_Evaluation';
    return `/api/grn_evaluation/${dataset.id}/${evalName}/image/${imageName}`;
  }

  /**
   * Register on-demand run to dataset
   */
  public registerOnDemandGrnRun(compare: boolean = false): void {
    const dataset = compare ? this.selectedDatasetCompare : this.selectedDataset;
    const selectedRun = compare ? this.selectedOnDemandRunCompare : this.selectedOnDemandRun;
    if (!dataset || !selectedRun) return;

    const formData = new FormData();
    formData.append('grn_evaluation_name', selectedRun);

    this.http.post(
      `/api/grn_evaluation/${dataset.id}/register_on_demand`,
      formData
    ).subscribe({
      next: () => {
        alert(`GRN run '${selectedRun}' registered successfully`);
      },
      error: (err) => {
        alert('Failed to register: ' + (err?.error?.detail || err.message));
      }
    });
  }

  // =========================================================================
  // GRN Evaluation Methods
  // =========================================================================
}



// CellGeometry/CellProperties/CellFeature/PropertyGroup/PropertyGroupItem moved to
// hexagon-view/cell-feature.types.ts (imported at the top of this file) so
// HexagonViewComponent can use the same types without a circular import back to this file.

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

interface tfGraphData {
  nodes: { id: string; name: string }[];
  links: { source: string; target: string; weight: number }[];
}

