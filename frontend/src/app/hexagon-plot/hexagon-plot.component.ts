import { Component, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import * as d3 from 'd3';
import * as Plotly from 'plotly.js-dist-min';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { HttpClient } from '@angular/common/http';
import { SessionService } from '../session.service';

@Component({
  selector: 'app-hexagon-plot',
  imports: [CommonModule, FormsModule, FilterableTableComponent],
  standalone: true,
  templateUrl: './hexagon-plot.component.html',
  styleUrl: './hexagon-plot.component.scss',
})
export class HexagonPlotComponent implements OnInit {
  constructor(
    private http: HttpClient,
    private sessionService: SessionService,
  ) {}

  // Define to use Math functions in the html template
  public Math = Math;

  // GeoJson
  public dataPath = 'assets/hexagons.geojson';
  public dataSetTitle =
    this.dataPath.split('/').pop()?.replace('.geojson', '') || 'Hexagon Plot';

  // Map svg and g elements
  private svg!: d3.Selection<SVGSVGElement, unknown, HTMLElement, any>;
  private g!: d3.Selection<SVGGElement, unknown, HTMLElement, any>;

  public selectedCell: CellFeature | null = null;
  public selectedCluster: number | null = null;
  public colorByProperty = 'cell_type';
  public selectedGeneSetGenie3: string | null = null;
  public selectedGeneSetSponge: string | null = null;
  public selectedRegulatoryScore: string | null = null;
  public selectedCellRegulatoryScores: { [key: string]: number } | null = null;
  private previousGeneSetGenie3: string | null = null;
  private previousGeneSetSponge: string | null = null;

  public selectedInterval: number = 0;
  public features: CellFeature[] = []; // public so that filterable table can update it
  public meta: { [key: string]: any } = {};

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

  public genie3Network: genie3RegGraphConnection[] = [];
  public spongeNetwork: spongeRegGraphConnection[] = [];
  public geneSetsGenie3: { [regulator: string]: string[] } = {};
  public geneSetsSponge: { [regulator: string]: string[] } = {};

  public minGenie3Edges: number = 25;
  public minSpongeEdges: number = 25;
  public genie3WeightCutoff: number | null = null;
  public spongePValueCutoff: number | null = null;
  public isLoadingSponge: boolean = false;
  public isLoadingGenie3: boolean = false;

  public coOccurrenceData: number[] = [];
  public coOccurrenceColumns: string[] = [];
  public coOccurrenceThreshold: number = 0.5;
  public maxInterval: number = 49;
  public clusterCount: number = 10;

  public colorableProperties = [
    'cell_type',
    'leiden',
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];
  public leidenCentralityProps = [
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
  ];
  public ligandReceptorScores: {
    [col: string]: { [index: string]: string | number };
  } | null = null;
  public showGlobalLigandReceptorScores: boolean = true;
  public showMoranI: boolean = true;
  public colorScale = d3
    .scaleOrdinal<string>()
    .range([
      '#FF7373',
      '#66cdaa',
      '#088da5',
      '#F0E442',
      '#0072B2',
      '#ffc3a0',
      '#CC79A7',
      '#E15759',
    ]);
  private continuousColorScale = d3.scaleSequential(d3.interpolateBlues);
  public currentLegendDomain: any[] = [];
  public currentLegendType: 'continuous' | 'categorical' = 'categorical';

  ngOnInit(): void {
    this.createHexagonPlot();
    this.loadAndRenderData(this.dataPath);
    this.updateHexColors();
  }

  private createHexagonPlot(): void {
    const width = 500;
    const height = 400;
    this.svg = d3
      .select('#hexbin')
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', [0, 0, 1000, 1000] as [number, number, number, number])
      .style('background-color', 'white')
      .style('overflow', 'hidden')
      .call(
        d3
          .zoom<SVGSVGElement, unknown>()
          .scaleExtent([1, 5])
          .extent([
            [0, 0],
            [width, height],
          ])
          .on('zoom', (event) => this.zoomed(event)),
      );

    this.g = this.svg.append('g');
  }

  private zoomed(event: d3.D3ZoomEvent<SVGSVGElement, unknown>): void {
    this.g.attr('transform', event.transform.toString());
  }

  private loadAndRenderData(dataPath: string): void {
    d3.json<GeoJsonData>(dataPath)
      .then((data) => {
        console.log('Data loaded:', data);

        if (!data) {
          throw new Error('Failed to load GeoJSON data');
        }

        // This is for showing all properties for coloring
        this.features = data.features;
        this.meta = data.meta;
        this.selectedRegulatoryScore =
          this.meta['grn_score_names']?.[0] || null;
        this.geneSetsGenie3 = this.meta['genie_genesets'] || {};
        this.geneSetsSponge = this.meta['sponge_genesets'] || {};
        this.selectedGeneSetGenie3 =
          Object.keys(this.meta['genie_genesets'] || {})[0] || null;
        this.selectedGeneSetSponge =
          Object.keys(this.meta['sponge_genesets'] || {})[0] || null;

        const firstProps = this.features[0]?.properties || {};
        this.colorableProperties = Object.keys(firstProps).filter((k) => {
          const val = firstProps[k];
          return typeof val === 'string' || typeof val === 'number';
        });

        this.colorableProperties.push('regulatory_scores');

        // sort in alphabetical order

        this.colorableProperties.sort((a, b) => a.localeCompare(b));

        const width = 1200;
        const height = 1000;

        const projection = d3.geoIdentity().fitSize([width, height], {
          type: 'FeatureCollection',
          features: data.features,
        });

        this.features = data.features;

        //this.colorScale.domain([
        //  ...new Set(
        //    this.features.map((f) =>
        //      String(f.properties[this.colorByProperty]),
        //    ),
        //  ),
        //]);
        this.colorScale.domain([
          ...new Set(
            this.features.map((f: CellFeature) => f.properties.cell_type),
          ),
        ]);
        this.currentLegendDomain = this.colorScale.domain();

        // Create a geoPath generator with the projection
        const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

        // Draw the map inside the zoomable group
        this.g
          .style('cursor', 'pointer')
          .append('g')
          .selectAll('path')
          .data(data.features)
          .join('path')
          .attr('d', (d) => pathGenerator(d))
          .attr('fill', (d) => {
            const total = d.properties.cell_type || 0;
            return this.colorScale(total.toString());
          })
          .style('opacity', 0.8)
          .on('mouseover', (event, d) => this.mouseOver(event, d))
          .on('mouseleave', (event, d) => this.mouseLeave(event, d))
          .on('click', (event, d) => this.openSidenav(event, d));

        this.colorScale.domain([
          ...new Set(
            data.features.map((f: CellFeature) => f.properties.cell_type),
          ),
        ]);
        this.renderLegend();
      })
      .catch((error) => {
        console.error('Error loading or rendering data:', error);
      });

    /* d3.csv('assets/sponge_network_smaller.csv', d3.autoType)
      .then((rows) => {
        console.log('WEird')
        this.spongeNetwork = rows.map((row) => ({
          source: String((row as any)['geneA'] ?? ''),
          target: String((row as any)['geneB'] ?? ''),
          p_adjusted: Number((row as any)['p.adj'] ?? 0),
        }));
        console.log('Sponge network loaded:', this.spongeNetwork);
      })
      .catch((error) => {
        console.error('Error loading sponge network:', error);
      });

    // Read genie3 csv
    d3.csv('assets/genie_network_filt.csv', d3.autoType)
      .then((rows) => {
        // Each row should have source, target, weight columns
        this.genie3Network = rows.map((row) => ({
          source: String((row as any)['regulatoryGene'] ?? ''),
          target: String((row as any)['targetGene'] ?? ''),
          weight: Number((row as any)['weight'] ?? 0),
        }));
        console.log('Genie3 network loaded:', this.genie3Network);
      })
      .catch((error) => {
        console.error('Error loading genie3 network:', error);
      });
 */
  }

  //public updateHexColors(): void {
  //  this.resetClusterExtension();
  //
  //  if (this.selectedCell && this.selectedCluster) {
  //    this.selectedCluster = null;
  //    this.clusterCells = [];
  //    this.clusterCellTypes = [];
  //    this.clusterCentralityAvg = {
  //      degree_centrality: 0,
  //      average_clustering: 0,
  //      closeness_centrality: 0,
  //    };
  //  }
  //
  //  if (this.selectedCell) {
  //    this.selectedCell = null;
  //  }
  //
  //  if (this.leidenCentralityProps.includes(this.colorByProperty)) {
  //    // Get all values for the selected centrality property
  //    const values = this.features.map(
  //      (f) => f.properties.leiden_centrality[this.colorByProperty],
  //    );
  //    const min = Math.min(...values);
  //    const max = Math.max(...values);
  //
  //    this.continuousColorScale.domain([min, max]);
  //
  //    // Update hexagon colors using the continuous scale
  //    this.g
  //      .selectAll<SVGPathElement, CellFeature>('path')
  //      .transition()
  //      .duration(300)
  //      .attr('fill', (d: CellFeature) =>
  //        this.continuousColorScale(
  //          d.properties.leiden_centrality[this.colorByProperty],
  //        ),
  //      );
  //  } else {
  //    // Categorical color scale for other properties
  //    this.colorScale.domain([
  //      ...new Set(
  //        this.features.map((f) => String(f.properties[this.colorByProperty])),
  //      ),
  //    ]);
  //    this.g
  //      .selectAll<SVGPathElement, CellFeature>('path')
  //      .transition()
  //      .duration(300)
  //      .style('stroke', 'transparent')
  //      .attr('fill', (d: CellFeature) =>
  //        this.colorScale(String(d.properties[this.colorByProperty])),
  //      );
  //  }
  //
  //  this.renderLegend();
  //}
  //

  private toNumber(v: unknown): number {
    if (typeof v === 'number') return v;
    if (typeof v === 'string') {
      const n = parseFloat(v);
      return Number.isFinite(n) ? n : NaN;
    }
    return NaN;
  }

  public onColorbyPropertyChange(): void {
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
        this.updateAucellGraphGenie3();
        this.updateAucellGraphSponge();
      } else if (
        this.selectedRegulatoryScore?.endsWith('sponge') &&
        this.selectedGeneSetSponge
      ) {
        this.fetchAndUpdate(
          this.selectedRegulatoryScore,
          this.selectedGeneSetSponge,
        );
        this.updateAucellGraphGenie3();
        this.updateAucellGraphSponge();
      }
    }
    this.updateHexColors();
  }

  public updateHexColors(): void {
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

    // 1) collect values (supports nested props like ligand_receptor_relationships)
    const valuesRaw = this.features.map((f) => {
      if (this.leidenCentralityProps.includes(this.colorByProperty)) {
        return f.properties.leiden_centrality[this.colorByProperty];
      }
      return f.properties[this.colorByProperty];
    });

    const numericValues = valuesRaw.map((v) => this.toNumber(v));
    const allNumbers = numericValues.every((n) => Number.isFinite(n));

    // Check if all values are integers (for categorical treatment)
    const allIntegers =
      allNumbers && numericValues.every((n) => Number.isInteger(n));

    // Check if we have a reasonable number of unique integer values for categorical treatment (here 20)
    const uniqueIntegerCount = allIntegers ? new Set(numericValues).size : 0;
    const shouldTreatAsCategorical = allIntegers && uniqueIntegerCount <= 20;

    const sel = this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .data(this.features);

    if (allNumbers && !shouldTreatAsCategorical && numericValues.length > 0) {
      // continuous scale - only if not integers or too many unique integers
      let min = Math.min(...numericValues);
      let max = Math.max(...numericValues);
      if (min === max) {
        const eps = min === 0 ? 1 : Math.abs(min) * 0.01;
        min -= eps;
        max += eps;
      }
      this.continuousColorScale.domain([min, max]);
      this.currentLegendDomain = [min, max];
      this.currentLegendType = 'continuous';

      sel
        .transition()
        .duration(300)
        .attr('stroke-width', 1)
        .attr('stroke', 'transparent')
        .attr('fill', (d) => {
          const raw = this.leidenCentralityProps.includes(this.colorByProperty)
            ? d.properties.leiden_centrality[this.colorByProperty]
            : d.properties[this.colorByProperty];
          const n = this.toNumber(raw);
          return Number.isFinite(n) ? this.continuousColorScale(n) : '#ccc';
        });
    } else {
      // categorical scale - for non-numeric, integers with few unique values, or mixed data
      const domain = [...new Set(valuesRaw.map((v: any) => String(v)))];
      this.colorScale.domain(domain);
      this.currentLegendDomain = domain;
      this.currentLegendType = 'categorical';

      sel
        .transition()
        .duration(300)
        .attr('stroke-width', 1)
        .attr('stroke', 'transparent')
        .attr('fill', (d) => {
          const raw = this.leidenCentralityProps.includes(this.colorByProperty)
            ? d.properties.leiden_centrality[this.colorByProperty]
            : d.properties[this.colorByProperty];
          return this.colorScale(String(raw));
        });
    }

    this.renderLegend();
  }

  public updateAucellGraphGenie3(): void {
    console.log('Updating AUCELL graph for Genie3...');
    this.isLoadingGenie3 = true;
    d3.select('#aucell_graph_genie3').selectAll('*').remove();

    if (!this.selectedGeneSetGenie3 || !this.genie3Network) {
      this.isLoadingGenie3 = false;
      return;
    }

    const regulator = this.selectedGeneSetGenie3;
    const targets = this.geneSetsGenie3[regulator] || [];

    const nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    const edges: { source: string; target: string; weight: number }[] = [];

    let candidateEdges: { source: string; target: string; weight: number }[] =
      [];

    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/geneset_connections_genie?gene_set_name=${encodeURIComponent(this.selectedGeneSetGenie3 ? this.selectedGeneSetGenie3 : '')}`,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          const data = res as {
            regulatoryGene: string;
            targetGene: string;
            weight: number;
          }[];
          console.log('Data', data);
          this.genie3Network = data.map((d) => ({
            source: d.regulatoryGene,
            target: d.targetGene,
            weight: d.weight,
          }));

          this.genie3Network.forEach((connection) => {
            if (
              connection.source === regulator ||
              connection.target === regulator ||
              targets.includes(connection.source) ||
              targets.includes(connection.target)
            ) {
              candidateEdges.push({
                source: connection.source,
                target: connection.target,
                weight: connection.weight,
              });
            }
          });

          console.log('Candidate edges before filtering:', candidateEdges);

          candidateEdges.sort((a, b) => b.weight - a.weight);
          candidateEdges = candidateEdges.slice(0, this.minGenie3Edges);

          this.genie3WeightCutoff =
            candidateEdges.length > 0
              ? Math.min(...candidateEdges.map((edge) => edge.weight))
              : null;

          // Create nodes from all edges (source and target)
          const nodeSet = new Set<string>();
          candidateEdges.forEach((edge) => {
            nodeSet.add(edge.source);
            nodeSet.add(edge.target);
          });
          // Add regulator to nodeset
          nodeSet.add(regulator);

          // For the nodes with top edges we reinset their original edges

          this.genie3Network.forEach((connection) => {
            if (
              nodeSet.has(connection.source) &&
              nodeSet.has(connection.target)
            ) {
              candidateEdges.push({
                source: connection.source,
                target: connection.target,
                weight: connection.weight,
              });
            }
          });

          // Create nodes array with proper groups
          Array.from(nodeSet).forEach((nodeId) => {
            if (nodeId === regulator) {
              if (!nodes.some((n) => n.id === nodeId)) {
                nodes.push({ id: nodeId, group: 0 }); // regulator
              }
            } else if (targets.includes(nodeId)) {
              if (!nodes.some((n) => n.id === nodeId)) {
                nodes.push({ id: nodeId, group: 1 }); // targets
              }
            } else {
              if (!nodes.some((n) => n.id === nodeId)) {
                nodes.push({ id: nodeId, group: 2 }); // neighbors
              }
            }
          });

          edges.push(...candidateEdges);

          // Create the graph
          const graph = {
            nodes: nodes.filter((node) => node.id && node.id.length > 0),
            edges: edges.filter(
              (edge) =>
                nodes.some((node) => node.id === edge.source) &&
                nodes.some((node) => node.id === edge.target),
            ),
          };

          // Create the graph visualization
          const width = 500;
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

          console.log('Network visualization complete');
          console.log('Genie3 Network:', this.genie3Network);
          console.log(
            `[Backend] Loaded Genie Connections for["${this.selectedGeneSetGenie3}]`,
          );
        },
        error: (err) =>
          console.error(
            `[Backend] Failed to load Genie Connections for["${this.selectedGeneSetGenie3}]`,
            err,
          ),
      });

    console.log('Genie3 Network 2:', this.genie3Network);
  }

  public updateAucellGraphSponge(): void {
    console.log('Updating AUCELL graph for Sponge...');
    d3.select('#aucell_graph_sponge').selectAll('*').remove();

    if (!this.selectedGeneSetSponge || !this.spongeNetwork) {
      return;
    }

    this.isLoadingSponge = true;
    const regulator = this.selectedGeneSetSponge;
    const targets = this.geneSetsSponge[regulator] || [];

    const nodes: { id: string; x?: number; y?: number; group: number }[] = [];
    const edges: { source: string; target: string; p_adjusted: number }[] = [];

    let candidateEdges: {
      source: string;
      target: string;
      p_adjusted: number;
    }[] = [];
    /*
    this.spongeNetwork.forEach(connection => {
     if ((connection.source === regulator || connection.target === regulator) ||
         targets.includes(connection.source) || targets.includes(connection.target)) {
       candidateEdges.push({
         source: connection.source,
         target: connection.target,
         p_adjusted: connection.p_adjusted
       });
     }
   }) */ this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/geneset_connections_sponge?gene_set_name=${encodeURIComponent(this.selectedGeneSetSponge ? this.selectedGeneSetSponge : '')}`,
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => {
          const data = res as {
            geneA: string;
            geneB: string;
            'p.adj': number;
            mscor: number;
          }[];

          this.isLoadingSponge = true;

          console.log('Sponge Network:', data);

          this.spongeNetwork = data.map((d) => ({
            source: d.geneA,
            target: d.geneB,
            p_adjusted: d['p.adj'],
            mscore: d['mscor'],
          }));

          // Push all edges of filtered Network (filtered to geneset edges)
          candidateEdges.push(
            ...this.spongeNetwork.filter((connection) => {
              return (
                connection.source === regulator ||
                connection.target === regulator ||
                targets.includes(connection.source) ||
                targets.includes(connection.target)
              );
            }),
          );
          // Sort such that smallest p_values are kept
          candidateEdges.sort((a, b) => a.p_adjusted - b.p_adjusted);
          candidateEdges = candidateEdges.slice(0, this.minSpongeEdges);

          this.spongePValueCutoff =
            candidateEdges.length > 0
              ? Math.max(...candidateEdges.map((edge) => edge.p_adjusted))
              : null;

          // Create nodes from all edges (source and target)
          const nodeSet = new Set<string>();
          candidateEdges.forEach((edge) => {
            nodeSet.add(edge.source);
            nodeSet.add(edge.target);
          });
          // Add regulator to nodeset
          nodeSet.add(regulator);

          this.spongeNetwork.forEach((connection) => {
            if (
              nodeSet.has(connection.source) &&
              nodeSet.has(connection.target)
            ) {
              candidateEdges.push({
                source: connection.source,
                target: connection.target,
                p_adjusted: connection.p_adjusted,
              });
            }
          });

          edges.push(...candidateEdges);

          // Create nodes array with proper groups
          Array.from(nodeSet).forEach((nodeId) => {
            if (nodeId === regulator) {
              if (!nodes.some((n) => n.id === nodeId)) {
                nodes.push({ id: nodeId, group: 0 }); // regulator
              }
            } else if (targets.includes(nodeId)) {
              if (!nodes.some((n) => n.id === nodeId)) {
                nodes.push({ id: nodeId, group: 1 }); // targets
              }
            } else {
              if (!nodes.some((n) => n.id === nodeId)) {
                nodes.push({ id: nodeId, group: 2 }); // neighbors
              }
            }
          });

          // Create the graph
          const graph = {
            nodes: nodes.filter((node) => node.id && node.id.length > 0),
            edges: edges.filter(
              (edge) =>
                nodes.some((node) => node.id === edge.source) &&
                nodes.some((node) => node.id === edge.target),
            ),
          };

          // Create the graph visualization
          const width = 500;
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
        },
        error: (err) =>
          console.error(
            `[Backend] Failed to load Sponge Connections for["${this.selectedGeneSetSponge}]`,
            err,
          ),
      });
  }

  public onMinEdgesChangeSponge(): void {
    if (this.selectedGeneSetSponge) {
      this.isLoadingSponge = true;
      setTimeout(() => this.updateAucellGraphSponge(), 50);
    }
  }

  public onMinEdgesChangeGenie3(): void {
    if (this.selectedGeneSetGenie3) {
      this.isLoadingGenie3 = true;
      setTimeout(() => this.updateAucellGraphGenie3(), 50);
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

    if (this.colorByProperty === 'leiden') {
      this.openClusterSidenav(cell.properties.leiden);
      this.extendCluster(cell.properties.leiden);
    } else {
      d3.select(event.target as SVGElement)
        .transition()
        .attr('stroke', 'black');
    }

    setTimeout(() => this.renderNhoodHeatmap(), 0);

    setTimeout(() => this.updateAucellGraphGenie3(), 0);
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
      setTimeout(() => this.renderNhoodHeatmap(), 100);
      setTimeout(() => this.updateAucellGraphGenie3(), 100);
    }
  }

  public onGeneSetChange(): void {
    // Check if Genie3 gene set has actually changed
    if (this.selectedGeneSetGenie3 !== this.previousGeneSetGenie3) {
      this.previousGeneSetGenie3 = this.selectedGeneSetGenie3;

      if (this.selectedGeneSetGenie3) {
        d3.select('#aucell_graph_genie3').selectAll('*').remove();
        this.isLoadingGenie3 = true;
        setTimeout(() => {
          this.updateAucellGraphGenie3();
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
          this.updateAucellGraphSponge();
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
    setTimeout(() => this.renderNhoodHeatmap(), 0);
    setTimeout(() => this.updateAucellGraphGenie3(), 0);
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

    if (this.clusterCells.length > 0) {
      const firstCell = this.clusterCells[0];
      this.clusterCentralityAvg = {
        degree_centrality:
          firstCell.properties.leiden_centrality['degree_centrality'] || 0,
        average_clustering:
          firstCell.properties.leiden_centrality['average_clustering'] || 0,
        closeness_centrality:
          firstCell.properties.leiden_centrality['closeness_centrality'] || 0,
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

  private renderNhoodHeatmap(): void {
    if (!this.selectedCell?.properties.leiden_nhood_enrichment) return;

    const enrichment = this.selectedCell.properties.leiden_nhood_enrichment;
    const leiden = this.selectedCell.properties.leiden;
    const n = enrichment.length;
    const clusterLabels = Array.from({ length: n }, (_, i) => `Cluster ${i}`);

    const data: Partial<Plotly.PlotData>[] = [
      {
        z: [enrichment],
        x: clusterLabels,
        y: [leiden.toString()],
        type: 'heatmap',
        colorscale: 'Viridis',
      },
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
        title: { text: '' },
        automargin: true,
        showticklabels: false,
        tickfont: { size: 10 },
      },
    };

    const container = document.getElementById('cluster-nhood-heatmap');
    if (!container) {
      console.error('Container cluster-nhood-heatmap not found');
      return;
    }

    if (!container) return;
    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, { displayModeBar: false });
  }

  public closeSidenav(): void {
    this.selectedCell = null;
    this.updateHexColors();
  }

  public updateCoOccurrenceTable(): void {
    if (this.features.length === 0 || this.selectedCluster === null) {
      this.coOccurrenceData = [];
      return;
    }

    // Get co-occurrence data from the first cell
    const firstCell = this.features[0];
    if (!firstCell?.properties?.leiden_co_occurrence) {
      console.warn('No leiden_co_occurrence data found');
      this.coOccurrenceData = [];
      return;
    }

    const coOccurrenceMatrix = firstCell.properties.leiden_co_occurrence;

    if (!Array.isArray(coOccurrenceMatrix)) {
      console.error(
        'Co-occurrence matrix is not an array:',
        coOccurrenceMatrix,
      );
      this.coOccurrenceData = [];
      return;
    }

    this.coOccurrenceData = [];
    console.log(coOccurrenceMatrix[1][this.selectedInterval]);
    try {
      for (let j = 0; j < this.clusterCount; j++) {
        if (
          Array.isArray(coOccurrenceMatrix[j]) &&
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

    // Create a color scale from light to dark based on value
    const maxValue = Math.max(...this.coOccurrenceData.flat());
    const intensity = Math.min(value / maxValue, 1);

    // Use a blue color scale
    const blue = Math.floor(255 - intensity * 200);
    const green = Math.floor(255 - intensity * 150);
    return `rgb(${blue}, ${green}, 255)`;
  }

  public getIntervalStats(): { min: number; max: number; avg: number } {
    const allValues = this.coOccurrenceData.flat().filter((val) => val > 0);
    if (allValues.length === 0) return { min: 0, max: 0, avg: 0 };

    const min = Math.min(...allValues);
    const max = Math.max(...allValues);
    const avg = allValues.reduce((sum, val) => sum + val, 0) / allValues.length;

    return { min, max, avg: Math.round(avg * 100) / 100 };
  }

  async fetchAndUpdate(columnName: string, index: string) {
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

  private renderLegend(): void {
    // Remove any existing legend
    this.svg.selectAll('.svg-legend').remove();

    if (this.currentLegendType === 'continuous') {
      // Continuous legend
      const [min, max] = this.currentLegendDomain as number[];
      const legendX = -100;
      const legendY = 50;
      const width = 250;
      const height = 30;
      const fontSize = 24;
      const padding = 15;

      // Create gradient for continuous legend
      const defs = this.svg.select('defs').empty()
        ? this.svg.append('defs')
        : this.svg.select('defs');

      defs.select('#svg-legend-gradient').remove();

      const gradient = defs
        .append('linearGradient')
        .attr('id', 'svg-legend-gradient')
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
          .attr('stop-color', this.continuousColorScale(value));
      }

      const legendG = this.svg
        .append('g')
        .attr('class', 'svg-legend')
        .attr('transform', `translate(${legendX},${legendY})`);

      // Measure title text width for dynamic background
      const titleText = this.colorByProperty.replace(/_/g, ' ');
      const tempSvg = this.svg.append('g').style('opacity', 0);
      const titleWidth =
        tempSvg
          .append('text')
          .text(titleText)
          .style('font-size', `${fontSize}px`)
          .style('font-weight', 'bold')
          .node()
          ?.getBBox().width || 0;

      // Measure min value
      const minText = min.toFixed(2);
      const minWidth =
        tempSvg
          .append('text')
          .text(minText)
          .style('font-size', `${fontSize}px`)
          .node()
          ?.getBBox().width || 0;

      // Measure max value
      const maxText = max.toFixed(2);
      const maxWidth =
        tempSvg
          .append('text')
          .text(maxText)
          .style('font-size', `${fontSize}px`)
          .node()
          ?.getBBox().width || 0;

      tempSvg.remove();

      // Calculate required dimensions
      const textHeight = fontSize * 1.2; // Approximate text height
      const requiredWidth = Math.max(
        width,
        titleWidth,
        minWidth + maxWidth + 20,
      );
      const bgWidth = requiredWidth + padding * 2;
      const bgHeight = height + textHeight * 2 + padding * 3;

      // Background
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

      // Gradient rectangle
      legendG
        .append('rect')
        .attr('x', (bgWidth - width) / 2 - padding)
        .attr('y', 0)
        .attr('width', width)
        .attr('height', height)
        .style('fill', 'url(#svg-legend-gradient)')
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

      // Title
      legendG
        .append('text')
        .attr('x', bgWidth / 2 - padding)
        .attr('y', -5)
        .attr('text-anchor', 'middle')
        .style('font-size', `${fontSize}px`)
        .style('font-weight', 'bold')
        .style('fill', '#333')
        .text(titleText);
    } else {
      // Categorical legend
      const categories = this.currentLegendDomain as string[];
      const legendX = -100;
      const legendY = 10;
      const itemHeight = 40;
      const rectHeight = 20;
      const rectWidth = 30;
      const fontSize = 24;
      const titlePadding = 15;

      // Create temporary text elements to measure actual width
      const tempSvg = this.svg.append('g').style('opacity', 0);

      // Measure title text
      const titleText = this.colorByProperty.replace(/_/g, ' ');
      const titleWidth =
        tempSvg
          .append('text')
          .text(titleText)
          .style('font-size', `${fontSize}px`)
          .style('font-weight', 'bold')
          .node()
          ?.getBBox().width || 0;

      // Measure category text widths
      const textNodes = tempSvg
        .selectAll('text')
        .data(categories)
        .enter()
        .append('text')
        .text((d) => d)
        .style('font-size', `${fontSize}px`);

      const maxTextWidth = Math.max(
        ...textNodes
          .nodes()
          .map((node) => (node as SVGGraphicsElement).getBBox().width),
      );
      tempSvg.remove();

      const itemWidth = Math.max(200, maxTextWidth + 60, titleWidth + 40);

      const legendG = this.svg
        .append('g')
        .attr('class', 'svg-legend')
        .attr('transform', `translate(${legendX},${legendY})`);

      // Calculate title height and total background height
      const titleHeight = fontSize + titlePadding;
      const backgroundHeight =
        categories.length * itemHeight + 20 + titleHeight;
      const backgroundWidth = itemWidth + 20;

      // Background - positioned to include title space
      legendG
        .append('rect')
        .attr('x', -10)
        .attr('y', -10)
        .attr('width', backgroundWidth)
        .attr('height', backgroundHeight)
        .style('fill', 'rgba(255, 255, 255, 0.9)')
        .attr('stroke', '#ccc')
        .attr('stroke-width', 1)
        .attr('rx', 5);

      // Add title
      legendG
        .append('text')
        .attr('x', backgroundWidth / 2 - 10)
        .attr('y', titlePadding + fontSize / 2)
        .attr('dy', '0.35em')
        .attr('text-anchor', 'middle')
        .style('font-size', `${fontSize}px`)
        .style('font-weight', 'bold')
        .style('fill', '#333')
        .text(titleText);

      categories.forEach((cat, i) => {
        const yPosition = i * itemHeight + titleHeight;
        const legendItem = legendG
          .append('g')
          .attr('transform', `translate(0, ${yPosition})`);

        // Color rectangle - centered vertically within the item height
        const rectY = (itemHeight - rectHeight) / 2;
        legendItem
          .append('rect')
          .attr('y', rectY)
          .attr('width', rectWidth)
          .attr('height', rectHeight)
          .style('fill', this.colorScale(cat))
          .attr('stroke', '#333')
          .attr('stroke-width', 0.5)
          .attr('rx', 2);

        // Text - aligned with the center of the rectangle
        const textY = rectY + rectHeight / 2;
        legendItem
          .append('text')
          .attr('x', rectWidth + 10)
          .attr('y', textY)
          .attr('dy', '0.35em')
          .style('font-size', `${fontSize}px`)
          .style('fill', '#333')
          .text(cat);

        // Add tooltip for full text
        legendItem.append('title').text(cat);
      });
    }
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
  leiden_nhood_enrichment: number[];
  leiden: number;
  color: string;
  aucell_genie3: { [key: string]: number };
  aucell_sponge: { [key: string]: number };
  leiden_centrality: { [key: string]: number };
  leiden_co_occurrence: number[][];
  [key: string]:
    | string
    | number
    | number[]
    | []
    | undefined
    | { [key: string]: any };
}

interface CellFeature {
  type: 'Feature';
  geometry: CellGeometry;
  properties: CellProperties;
}

interface GeoJsonData {
  type: 'FeatureCollection';
  features: CellFeature[];
  meta: { [key: string]: any };
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
