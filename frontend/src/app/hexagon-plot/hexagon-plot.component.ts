import { Component, OnInit, ɵisComponentDefPendingResolution } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import * as d3 from 'd3';
import * as Plotly from 'plotly.js-dist-min';

@Component({
  selector: 'app-hexagon-plot',
  imports: [CommonModule, FormsModule],
  standalone: true,
  templateUrl: './hexagon-plot.component.html',
  styleUrl: './hexagon-plot.component.scss'
})
export class HexagonPlotComponent implements OnInit {
  // Define to use Math functions in the html template
  public Math = Math;
  // GeoJson
  public dataPath = 'assets/hexagons_GSM6592049_M2_aucell.geojson';
  public dataSetTitle = this.dataPath.split('/').pop()?.replace('.geojson', '') || 'Hexagon Plot';

  // Map svg and g elements
  private svg!: d3.Selection<SVGSVGElement, unknown, HTMLElement, any>;
  private g!: d3.Selection<SVGGElement, unknown, HTMLElement, any>;

  public selectedCell: CellFeature | null = null;
  public selectedCluster: number | null = null;
  public colorByProperty = 'cell_type';
  public selectedGeneSetGenie3: string | null = null;
  public selectedGeneSetSponge: string | null = null;
  private previousGeneSetGenie3: string | null = null;
  private previousGeneSetSponge: string | null = null;
  public selectedCellAssociatedGeneSetsGenie3: string[] = [];
  public selectedCellAssociatedGeneSetsSponge: string[] = [];
  public selectedInterval: number = 0;
  private features: CellFeature[] = [];

  public clusterCells: CellFeature[] = [];
  public clusterCellTypes: { type: string, count: number, percentage: string }[] = [];
  public clusterCentralityAvg: { degree_centrality: number, average_clustering: number, closeness_centrality: number } = {
    degree_centrality: 0,
    average_clustering: 0,
    closeness_centrality: 0
  };

  public genie3Network: genie3RegGraphConnection[] = [];
  public spongeNetwork: spongeRegGraphConnection[] = [];
  public geneSetsGenie3: { [regulator: string]: string[] } = {};
  public geneSetsSponge: { [regulator: string]: string[] } = {};

  public maxSpongeNodes: number = 25;
  public maxGenie3Nodes: number = 25;
  public spongePValueCutoff: number | null = null;
  public genie3WeightCutoff: number | null = null;
  public isLoadingSponge: boolean = false;
  public isLoadingGenie3: boolean = false;

  public coOccurrenceData: number[] = [];
  public coOccurrenceColumns: string[] = [];
  public coOccurrenceThreshold: number = 0.5;
  public maxInterval: number = 49;
  public clusterCount: number = 10;

  public colorableProperties = ['cell_type', 'leiden', 'degree_centrality', 'average_clustering', 'closeness_centrality'];
  public leidenCentralityProps = ["degree_centrality", "average_clustering", "closeness_centrality"];

  public colorScale = d3
    .scaleOrdinal<string>()
    .range([
      "#FF7373",
      "#66cdaa",
      "#088da5",
      "#F0E442",
      "#0072B2",
      "#ffc3a0",
      "#CC79A7",
      "#E15759"
    ])
  private continuousColorScale = d3.scaleSequential(d3.interpolateBlues);


  constructor() {

  }
  ngOnInit(): void {
    this.createHexagonPlot();
    this.loadAndRenderData(this.dataPath);

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
          .extent([[0, 0], [width, height]])
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

        const width = 1200;
        const height = 1000;


        const projection = d3.geoIdentity().fitSize([width, height], {
          type: 'FeatureCollection',
          features: data.features,
        });

        this.features = data.features;

        this.colorScale.domain([...new Set(this.features.map(f => String(f.properties[this.colorByProperty])))]);


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
          .on('click', (event, d) => this.openSidenav(event, d))

        this.colorScale.domain([...new Set(data.features.map((f: CellFeature) => f.properties.cell_type))]);
        this.renderLegend();
      })
      .catch((error) => {
        console.error('Error loading or rendering data:', error);
      });

    // Read genie3 csv
    d3.csv('assets/genie3_BRCA_mrn.top_100k.csv', d3.autoType).then((rows) => {
      // Each row should have source, target, weight columns
      this.genie3Network = rows.map(row => ({
        source: String((row as any)['regulatoryGene'] ?? ''),
        target: String((row as any)['targetGene'] ?? ''),
        weight: Number((row as any)['weight'] ?? 0)
      }));
      console.log('Genie3 network loaded:', this.genie3Network);
    }).catch((error) => {
      console.error('Error loading genie3 network:', error);
    });

    // Read sponge csv filtered
    d3.tsv('assets/breast_invasive_carcinoma_filtered_Network.tsv', d3.autoType).then((rows) => {
      // Each row should have source, target, p_adjusted columns
      this.spongeNetwork = rows.map(row => ({
        source: String((row as any)['geneA'] ?? ''),
        target: String((row as any)['geneB'] ?? ''),
        p_adjusted: Number((row as any)['p.adj'] ?? 0)
      }));
      console.log('Sponge network loaded:', this.spongeNetwork);
    }).catch((error) => {
      console.error('Error loading sponge network:', error);
    });

    // Read genie3 gene sets
    d3.json<{ [regulator: string]: string[] }>('assets/genie3_gene_sets_GSM6592049_M2.json')
      .then((data) => {
        this.geneSetsGenie3 = data || {};
        console.log('Genie3 gene sets loaded:', this.geneSetsGenie3);
      })
      .catch((error) => {
        console.error('Error loading genie3 gene sets:', error);
      });

    // Read sponge gene sets
    d3.json<{ [regulator: string]: string[] }>('assets/sponge_gene_sets_GSM6592049_M2.json')
      .then((data) => {
        this.geneSetsSponge = data || {};
        console.log('Sponge gene sets loaded:', this.geneSetsSponge);
      })
      .catch((error) => {
        console.error('Error loading sponge gene sets:', error);
      });

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
        closeness_centrality: 0
      };

    }

    if (this.selectedCell) {
      this.selectedCell = null;
    }


    if (this.leidenCentralityProps.includes(this.colorByProperty)) {
      // Get all values for the selected centrality property
      const values = this.features.map(f => f.properties.leiden_centrality[this.colorByProperty]);
      const min = Math.min(...values);
      const max = Math.max(...values);

      this.continuousColorScale.domain([min, max]);

      // Update hexagon colors using the continuous scale
      this.g.selectAll<SVGPathElement, CellFeature>('path')
        .transition()
        .duration(300)
        .attr('fill', (d: CellFeature) =>
          this.continuousColorScale(d.properties.leiden_centrality[this.colorByProperty])
        );
    } else {
      // Categorical color scale for other properties
      this.colorScale.domain([...new Set(this.features.map(f => String(f.properties[this.colorByProperty])))]);
      this.g.selectAll<SVGPathElement, CellFeature>('path')
        .transition()
        .duration(300)
        .style('stroke', 'transparent')
        .attr('fill', (d: CellFeature) => this.colorScale(String(d.properties[this.colorByProperty])));
    }

    this.renderLegend();
  }

  public onMaxNodesChangeSponge(): void {
    console.log('Max nodes changed to:', this.maxSpongeNodes);

    // Only update if a sponge gene set is selected
    if (this.selectedGeneSetSponge) {
      this.isLoadingSponge = true;
      setTimeout(() => this.updateAucellGraphSponge(), 50);
    }
  }

  public onMaxNodesChangeGenie3(): void {
    console.log('Max nodes changed to:', this.maxGenie3Nodes);
    // Only update if a Genie3 gene set is selected
    if (this.selectedGeneSetGenie3) {
      this.isLoadingGenie3 = true;
      setTimeout(() => this.updateAucellGraphGenie3(), 50);
    }
  }

  public onGeneSetChange(): void {
    console.log('Genie3 Gene set:', this.selectedGeneSetGenie3);
    console.log('Sponge Gene set:', this.selectedGeneSetSponge);

    // Check if Genie3 gene set has actually changed
    if (this.selectedGeneSetGenie3 !== this.previousGeneSetGenie3) {
      this.previousGeneSetGenie3 = this.selectedGeneSetGenie3;

      if (this.selectedGeneSetGenie3) {
        d3.select('#aucell_graph_genie3').selectAll('*').remove();
        this.isLoadingGenie3 = true;
        setTimeout(() => {
          this.updateAucellGraphGenie3();
        }, 100);
      } else {
        // Clear Genie3 graph if no gene set is selected
        d3.select('#aucell_graph_genie3').selectAll('*').remove();
        this.isLoadingGenie3 = false;
      }
    }

    // Check if Sponge gene set has actually changed
    if (this.selectedGeneSetSponge !== this.previousGeneSetSponge) {
      this.previousGeneSetSponge = this.selectedGeneSetSponge;

      if (this.selectedGeneSetSponge) {
        d3.select('#aucell_graph_sponge').selectAll('*').remove();
        console.log('Updating Sponge graph for:', this.selectedGeneSetSponge);
        console.log('Sponge targets available:', this.geneSetsSponge[this.selectedGeneSetSponge]?.length || 0);
        this.isLoadingSponge = true;
        setTimeout(() => {
          this.updateAucellGraphSponge();
        }, 100);
      } else {
        // Clear Sponge graph if no gene set is selected
        d3.select('#aucell_graph_sponge').selectAll('*').remove();
        this.isLoadingSponge = false;
      }
    }
  }

  public updateAucellGraphGenie3(): void {

    d3.select('#aucell_graph_genie3').selectAll('*').remove();


    if (!this.selectedGeneSetGenie3 || !this.genie3Network) {
      this.isLoadingGenie3 = false;
      return;
    }

    const regulator = this.selectedGeneSetGenie3;
    const targets = (this.geneSetsGenie3)[regulator] || [];


    const nodes: { id: string, x?: number, y?: number, group: number }[] = [];
    const edges: { source: string, target: string, weight: number }[] = [];

    let candidateEdges: { source: string, target: string, weight: number }[] = [];

    // Helper function to get weight from Genie3 network
    const getGenie3Weight = (source: string, target: string): number => {
      const connection = this.genie3Network.find(conn =>
        (conn.source === source && conn.target === target) ||
        (conn.source === target && conn.target === source)
      );
      return connection ? connection.weight : 0;
    };

    // Regulator node is always added
    nodes.push({ id: regulator, group: 0 });

    // Add target edges with actual weights
    targets.forEach(target => {
      const weight = getGenie3Weight(regulator, target);
      candidateEdges.push({ source: regulator, target, weight });
    });

    // If network is small, find neighbors and add their edges to candidateEdges
    if (candidateEdges.length < this.maxGenie3Nodes) {
      const neighborSet = new Set<string>();
      const allMainNodes = [regulator, ...targets];

      // Find all neighbors connected to main nodes
      allMainNodes.forEach((mainNode: string) => {
        this.genie3Network.forEach(connection => {
          if (connection.source === mainNode && !allMainNodes.includes(connection.target)) {
            neighborSet.add(connection.target);
          } else if (connection.target === mainNode && !allMainNodes.includes(connection.source)) {
            neighborSet.add(connection.source);
          }
        });
      });

      // Add all actual edges involving neighbors
      const allPotentialNodes = [...allMainNodes, ...Array.from(neighborSet)];

      this.genie3Network.forEach(connection => {
        if (allPotentialNodes.includes(connection.source) && allPotentialNodes.includes(connection.target)) {
          // Check if this edge is not already in candidateEdges
          const edgeExists = candidateEdges.some(edge =>
            (edge.source === connection.source && edge.target === connection.target) ||
            (edge.source === connection.target && edge.target === connection.source)
          );

          if (!edgeExists) {
            candidateEdges.push({
              source: connection.source,
              target: connection.target,
              weight: connection.weight
            });
          }
        }
      });
    }

    // Sort candidateEdges by weight in descending order (highest weights first)
    candidateEdges.sort((a, b) => b.weight - a.weight);
    // Only keep the top N edges by weight
    candidateEdges = candidateEdges.slice(0, this.maxGenie3Nodes);

    // Create nodes from all edges (source and target)
    const nodeSet = new Set<string>();
    candidateEdges.forEach(edge => {
      nodeSet.add(edge.source);
      nodeSet.add(edge.target);
    });

    // Ensure regulator is always included
    nodeSet.add(regulator);

    // Create nodes array with proper groups
    Array.from(nodeSet).forEach(nodeId => {
      if (nodeId === regulator) {
        if (!nodes.some(n => n.id === nodeId)) {
          nodes.push({ id: nodeId, group: 0 }); // regulator
        }
      } else if (targets.includes(nodeId)) {
        if (!nodes.some(n => n.id === nodeId)) {
          nodes.push({ id: nodeId, group: 1 }); // targets
        }
      } else {
        if (!nodes.some(n => n.id === nodeId)) {
          nodes.push({ id: nodeId, group: 2 }); // neighbors
        }
      }
    });

    // Add all candidate edges to final edges
    edges.push(...candidateEdges);



    // Create the graph object
    const graph = {
      nodes: nodes.filter(node => node.id && node.id.length > 0),
      edges: edges.filter(edge =>
        nodes.some(node => node.id === edge.source) &&
        nodes.some(node => node.id === edge.target)
      )
    };

    // Get Genie3 weight cutoff: For current edges in the graph what is the max weight (excluding self-loops)
    const nonSelfLoopEdges = edges.filter(edge => edge.source !== edge.target);
    this.genie3WeightCutoff = nonSelfLoopEdges.length > 0 ? Math.max(...nonSelfLoopEdges.map(edge => edge.weight)) : null;

    // Create the graph visualization
    const width = 500;
    const height = 300;
    const textPadding = 40;

    const svg = d3.select('#aucell_graph_genie3')
      .append('svg')
      .attr('width', width)
      .attr('height', height + textPadding)
      .style('background-color', '#f8f9fa');

    // Draw links (edges)
    const link = svg.append('g')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.6)
      .selectAll('line')
      .data(graph.edges)
      .enter()
      .append('line')
      .attr('stroke-width', (d: any) => Math.max(1, Math.sqrt(d.weight) * 10))
      .attr('stroke', (d: any) => {
        // Color edges based on weight
        const intensity = Math.min(d.weight * 10, 1);
        return d3.interpolateReds(0.3 + intensity * 0.7);
      });

    // Draw nodes
    const node = svg.append('g')
      .attr('stroke', '#fff')
      .attr('stroke-width', 1.5)
      .selectAll('circle')
      .data(graph.nodes)
      .enter()
      .append('circle')
      .attr('r', (d: any) => {
        switch (d.group) {
          case 0: return 15;
          case 1: return 12;
          case 2: return 8;
          case 3: return 6;
          default: return 10;
        }
      })
      .attr('fill', (d: any) => {
        switch (d.group) {
          case 0: return '#e41a1c';
          case 1: return '#377eb8';
          case 2: return '#4daf4a';
          case 3: return '#ff7f00';
          default: return '#999';
        }
      })

    // Add labels
    const labels = svg.append('g')
      .selectAll('text')
      .data(graph.nodes)
      .enter()
      .append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '.35em')
      .style('font-size', (d: any) => d.group === 0 ? '12px' : '10px')
      .style('font-weight', (d: any) => d.group === 0 ? 'bold' : 'normal')
      .style('fill', '#333')
      .text((d: any) => d.id.length > 8 ? d.id.substring(0, 8) + '...' : d.id);

    const simulation = d3.forceSimulation(graph.nodes)
      .force('link', d3.forceLink(graph.edges).id((d: any) => d.id).distance(10).strength(0.5))
      .force('charge', d3.forceManyBody().strength(-500))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(15))
      .force('boundary', () => {
        graph.nodes.forEach((node: any) => {
          const radius = node.group === 0 ? 15 : node.group === 1 ? 12 : node.group === 2 ? 8 : 6;

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

      node
        .attr('cx', (d: any) => d.x)
        .attr('cy', (d: any) => d.y);

      labels
        .attr('x', (d: any) => d.x)
        .attr('y', (d: any) => d.y);

    });

    setTimeout(() => {
      simulation.stop();
      this.isLoadingGenie3 = false;
    }, 2000);

    // AUCell Score Label
    svg.append('text')
      .attr('x', width / 2)
      .attr('y', height + textPadding - 40)
      .attr('text-anchor', 'middle')
      .style('font-size', '14px')
      .style('fill', '#333')
      .text(`AUCell Score ${this.selectedGeneSetGenie3}: ${this.selectedCell?.properties.aucell_genie3?.[this.selectedGeneSetGenie3]?.toFixed(3) || 'N/A'}`);

    svg.append('text')
      .attr('x', width / 2)
      .attr('y', height + textPadding - 25)
      .attr('text-anchor', 'middle')
      .style('font-size', '12px')
      .style('fill', '#666')
      .text(`Edges: ${edges.length}, Nodes: ${nodes.length}, Max Weight: ${this.genie3WeightCutoff ? this.genie3WeightCutoff.toFixed(3) : 'N/A'}`);

  }

  private handleNodeClick(nodeData: any): void {
    console.log('Node clicked:', nodeData.id);

    // Get all node IDs from the current graph
    const currentGraphNodes = this.getCurrentGraphNodes();

    // Generate gProfiler URL for the current set of nodes
    const gProfilerUrl = this.generateGProfilerUrl(currentGraphNodes);

    if (gProfilerUrl) {
      window.open(gProfilerUrl, '_blank');
    } else {
      console.warn('Could not generate gProfiler URL');
    }
  }

  private getCurrentGraphNodes(): string[] {
    if (!this.selectedGeneSetGenie3 || !this.geneSetsGenie3) {
      return [];
    }

    const regulator = this.selectedGeneSetGenie3;
    const targets = this.geneSetsGenie3[regulator] || [];


    return [regulator, ...targets];
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
      console.warn('Could not generate gProfiler URL for Genie3 gene set:', this.selectedGeneSetGenie3);
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


  public updateAucellGraphSponge(): void {

    // Reset p-value cutoff
    this.spongePValueCutoff = null;

    // Clear previous graph
    d3.select('#aucell_graph_sponge').selectAll('*').remove();

    if (!this.selectedGeneSetSponge || !this.spongeNetwork || this.spongeNetwork.length === 0) {
      console.log('Missing selectedGeneSetSponge or spongeNetwork');
      this.isLoadingSponge = false;
      return;
    }

    try {
      const regulator = this.selectedGeneSetSponge;
      const targets = (this.geneSetsSponge)[regulator] || [];

      console.log('Sponge Regulator:', regulator);
      console.log('Sponge Targets:', targets);
      console.log('Max nodes allowed:', this.maxSpongeNodes);

      // Use slider value for MAX_NODES
      const MAX_NODES = this.maxSpongeNodes;
      const MAX_EDGES = Math.min(100, MAX_NODES * 2); // Scale edges with nodes
      const pThreshold = 0.01;

      const nodes: { id: string, x?: number, y?: number, group: number }[] = [];
      let allPotentialEdges: { source: string, target: string, p_adjusted: number }[] = [];

      // Add regulator node
      nodes.push({ id: regulator, group: 0 });

      // Add target nodes (limited by slider)
      const maxTargets = Math.floor(MAX_NODES * 0.6); // 60% of max nodes for targets
      const limitedTargets = targets.slice(0, Math.min(targets.length, maxTargets));
      limitedTargets.forEach((target: string) => {
        nodes.push({ id: target, group: 1 });
      });

      // Pre-filter sponge network
      const relevantConnections = this.spongeNetwork.filter(connection =>
        connection.p_adjusted < pThreshold &&
        (connection.source === regulator || limitedTargets.includes(connection.source) ||
          connection.target === regulator || limitedTargets.includes(connection.target))
      );

      console.log(`Filtered connections from ${this.spongeNetwork.length} to ${relevantConnections.length}`);

      // Get direct regulator-target edges
      relevantConnections.forEach(connection => {
        if ((connection.source === regulator && limitedTargets.includes(connection.target)) ||
          (connection.target === regulator && limitedTargets.includes(connection.source))) {
          allPotentialEdges.push({
            source: connection.source,
            target: connection.target,
            p_adjusted: connection.p_adjusted
          });
        }
      });

      // Add neighbors if space available
      if (nodes.length < MAX_NODES) {
        const neighborSet = new Set<string>();
        const allMainNodes = [regulator, ...limitedTargets];

        const neighborConnections: Array<{ gene: string, p_adjusted: number }> = [];

        allMainNodes.forEach((mainNode: string) => {
          relevantConnections.forEach(connection => {
            if (connection.source === mainNode && !allMainNodes.includes(connection.target)) {
              neighborConnections.push({ gene: connection.target, p_adjusted: connection.p_adjusted });
            } else if (connection.target === mainNode && !allMainNodes.includes(connection.source)) {
              neighborConnections.push({ gene: connection.source, p_adjusted: connection.p_adjusted });
            }
          });
        });

        // Sort by significance and take only top neighbors
        neighborConnections.sort((a, b) => a.p_adjusted - b.p_adjusted);
        const maxNeighbors = Math.min(neighborConnections.length, MAX_NODES - nodes.length);

        neighborConnections.slice(0, maxNeighbors).forEach(({ gene }) => {
          if (!neighborSet.has(gene)) {
            neighborSet.add(gene);
            nodes.push({ id: gene, group: 2 });
          }
        });

        console.log('Added neighbors:', Array.from(neighborSet));

        // Add edges involving neighbors
        const allNodeIds = new Set(nodes.map(n => n.id));
        relevantConnections.forEach(connection => {
          const sourceInNetwork = allNodeIds.has(connection.source);
          const targetInNetwork = allNodeIds.has(connection.target);

          if (sourceInNetwork && targetInNetwork) {
            const edgeExists = allPotentialEdges.some(e =>
              (e.source === connection.source && e.target === connection.target) ||
              (e.source === connection.target && e.target === connection.source)
            );

            if (!edgeExists) {
              allPotentialEdges.push({
                source: connection.source,
                target: connection.target,
                p_adjusted: connection.p_adjusted
              });
            }
          }
        });
      }

      // Filter to most significant edges
      allPotentialEdges.sort((a, b) => a.p_adjusted - b.p_adjusted);
      const edgeLimit = Math.min(allPotentialEdges.length, MAX_EDGES);
      const edges = allPotentialEdges.slice(0, edgeLimit);

      // Set the p-value cutoff (the highest p-value that made it into the network)
      if (edges.length > 0) {
        this.spongePValueCutoff = edges[edges.length - 1].p_adjusted;
      }

      console.log(`Using top ${edges.length} edges out of ${allPotentialEdges.length}`);
      console.log('P-value cutoff:', this.spongePValueCutoff);

      // Remove disconnected nodes, but ALWAYS keep the regulator (group 0)
      const connectedNodeIds = new Set<string>();
      edges.forEach(edge => {
        connectedNodeIds.add(edge.source);
        connectedNodeIds.add(edge.target);
      });

      // Always include the regulator node, even if it has no connections
      connectedNodeIds.add(regulator);

      const filteredNodes = nodes.filter(node =>
        node.group === 0 || connectedNodeIds.has(node.id)
      );

      console.log('Final nodes:', filteredNodes.length, 'Final edges:', edges.length);

      // Since regulator is always kept, check if we have any meaningful network
      if (filteredNodes.length <= 1) {
        const svg = d3.select('#aucell_graph_sponge')
          .append('svg')
          .attr('width', 500)
          .attr('height', 300);

        svg.append('text')
          .attr('x', 250)
          .attr('y', 150)
          .attr('text-anchor', 'middle')
          .style('font-size', '14px')
          .style('fill', '#666')
          .text(`No significant connections for ${regulator}`);

        this.isLoadingSponge = false;
        return;
      }

      // Create the graph
      const graph = {
        nodes: filteredNodes,
        edges: edges
      };

      // Create visualization
      const width = 500;
      const height = 300;
      const textPadding = 60; // Increased for p-value info

      const svg = d3.select('#aucell_graph_sponge')
        .append('svg')
        .attr('width', width)
        .attr('height', height + textPadding)
        .style('background-color', '#f8f9fa');

      // Simplified links
      const link = svg.append('g')
        .attr('stroke', '#999')
        .attr('stroke-opacity', 0.6)
        .selectAll('line')
        .data(graph.edges)
        .enter()
        .append('line')
        .attr('stroke-width', 2)
        .attr('stroke', '#666');

      // Simplified nodes
      const node = svg.append('g')
        .attr('stroke', '#fff')
        .attr('stroke-width', 2)
        .selectAll('circle')
        .data(graph.nodes)
        .enter()
        .append('circle')
        .attr('r', (d: any) => {
          switch (d.group) {
            case 0: return 12; // regulator
            case 1: return 10; // targets
            case 2: return 8;  // neighbors
            default: return 8;
          }
        })
        .attr('fill', (d: any) => {
          switch (d.group) {
            case 0: return '#e41a1c'; // regulator - red
            case 1: return '#377eb8'; // targets - blue
            case 2: return '#4daf4a'; // neighbors - green
            default: return '#999';
          }
        })
        .style('cursor', 'pointer')
        .on('click', (event: MouseEvent, d: any) => {
          this.handleNodeClick(d);
        })
        .on('mouseover', function (event: MouseEvent, d: any) {
          d3.select(this).attr('stroke', '#333').attr('stroke-width', 3);
        })
        .on('mouseout', function (event: MouseEvent, d: any) {
          d3.select(this).attr('stroke', '#fff').attr('stroke-width', 2);
        });

      // Simplified labels
      const labels = svg.append('g')
        .selectAll('text')
        .data(graph.nodes)
        .enter()
        .append('text')
        .attr('text-anchor', 'middle')
        .attr('dy', '.35em')
        .style('font-size', '10px')
        .style('fill', '#333')
        .style('pointer-events', 'none')
        .text((d: any) => d.id.length > 6 ? d.id.substring(0, 6) + '...' : d.id);

      // Add tooltips
      node.append('title')
        .text((d: any) => `${d.id}\nGroup: ${d.group === 0 ? 'Regulator' : d.group === 1 ? 'Target' : 'Neighbor'}\nClick to analyze in gProfiler`);

      // Simulation
      const simulation = d3.forceSimulation(graph.nodes)
        .alphaDecay(0.05)
        .force('link', d3.forceLink(graph.edges).id((d: any) => d.id).distance(30).strength(0.3))
        .force('charge', d3.forceManyBody().strength(-200))
        .force('center', d3.forceCenter(width / 2, height / 2))
        .force('collision', d3.forceCollide().radius(15));

      // Throttled tick updates
      let tickCount = 0;
      simulation.on('tick', () => {
        tickCount++;
        if (tickCount % 3 === 0) {
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
        }
      });


      setTimeout(() => {
        simulation.stop();
        this.isLoadingSponge = false;
      }, 2000);

      // Add AUCell Score
      svg.append('text')
        .attr('x', width / 2)
        .attr('y', height + textPadding - 40)
        .attr('text-anchor', 'middle')
        .style('font-size', '14px')
        .style('fill', '#333')
        .text(`AUCell Score ${this.selectedGeneSetSponge}: ${this.selectedCell?.properties.aucell_sponge?.[this.selectedGeneSetSponge]?.toFixed(3) || 'N/A'}`);

      // Add network info
      svg.append('text')
        .attr('x', width / 2)
        .attr('y', height + textPadding - 25)
        .attr('text-anchor', 'middle')
        .style('font-size', '12px')
        .style('fill', '#666')
        .text(`${edges.length} connections, ${filteredNodes.length} nodes`);

      console.log('Optimized sponge network visualization complete');
    } catch (error) {
      console.error('Error generating Sponge graph:', error);
      this.isLoadingSponge = false;
    }
  }

  private mouseOver(event: MouseEvent, d: CellFeature): void {
    d3.selectAll('.Country')
      .transition()
      .duration(200)
      .style('opacity', 0.5)
      .style('stroke', 'transparent');

    d3.select(event.target as SVGElement)
      .transition()
      .duration(200)
      .style('opacity', 0.8)
      .style('stroke', 'black');
  }

  private getAssociatedGeneSets(cell: CellFeature): string[][] {
    const genie3GeneSets = Object.keys(cell.properties.aucell_genie3 || {})
      .filter(key => cell.properties.aucell_genie3[key] > 0.6);
    const spongeGeneSets = Object.keys(cell.properties.aucell_sponge || {})
      .filter(key => cell.properties.aucell_sponge[key] > 0.6);
    console.log('Associated Genie3 gene sets:', genie3GeneSets);
    console.log('Associated Sponge gene sets:', spongeGeneSets);

    return [genie3GeneSets, spongeGeneSets];
  }


  private mouseLeave(event: MouseEvent, d: CellFeature): void {
    if (this.selectedCell && d.properties.barcode === this.selectedCell.properties.barcode) return;
    d3.selectAll('.Country')
      .transition()
      .duration(200)
      .style('opacity', 0.8)
      .style('stroke', 'transparent');

    d3.select(event.target as SVGElement)
      .transition()
      .duration(200)
      .style('stroke', 'transparent');
  }

  // Regular Cell Specific Sidenav

  public openSidenav(event: MouseEvent, cell: CellFeature): void {

    this.selectedCell = cell;
    console.log('Selected cell:', this.selectedCell);

    [this.selectedCellAssociatedGeneSetsGenie3, this.selectedCellAssociatedGeneSetsSponge] = this.getAssociatedGeneSets(cell);

    this.selectedGeneSetGenie3 = null;
    this.selectedGeneSetSponge = null;
    this.previousGeneSetGenie3 = null;
    this.previousGeneSetSponge = null;

    if (this.colorByProperty === 'leiden') {
      this.openClusterSidenav(cell.properties.leiden);
      this.extendCluster(cell.properties.leiden);
      setTimeout(() => this.renderNhoodHeatmap(), 0);
    }
    else {
      // Remove Black outline from any previously selected cell
      this.g.selectAll<SVGPathElement, CellFeature>('path')
        .style('stroke', 'transparent');

      // Highlight the selected cell
      d3.select(event.target as SVGElement)
        .transition()
        .style('stroke', 'black');
    }
  }

  // Clusterview Specific Sidenav

  public openClusterSidenav(clusterId: number): void {
    this.selectedCluster = clusterId;
    this.clusterCells = this.features.filter(cell => cell.properties.leiden === clusterId);
    this.calculateClusterStats();

    // Initialize co-occurrence table for this cluster
    this.updateCoOccurrenceTable();

    if (this.clusterCells.length > 0) {
      this.selectedCell = this.clusterCells[0];
      setTimeout(() => this.renderNhoodHeatmap(), 100);
    }

  }

  public selectCellFromCluster(cell: CellFeature): void {
    this.selectedCell = cell;
    setTimeout(() => this.renderNhoodHeatmap(), 0);
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
    this.clusterCells.forEach(cell => {
      const cellType = cell.properties.cell_type;
      cellTypeMap.set(cellType, (cellTypeMap.get(cellType) || 0) + 1);
    });

    this.clusterCellTypes = Array.from(cellTypeMap.entries())
      .map(([type, count]) => ({
        type,
        count,
        percentage: ((count / this.clusterCells.length) * 100).toFixed(1)
      }))
      .sort((a, b) => b.count - a.count);

    if (this.clusterCells.length > 0) {
      const firstCell = this.clusterCells[0];
      this.clusterCentralityAvg = {
        degree_centrality: firstCell.properties.leiden_centrality['degree_centrality'] || 0,
        average_clustering: firstCell.properties.leiden_centrality['average_clustering'] || 0,
        closeness_centrality: firstCell.properties.leiden_centrality['closeness_centrality'] || 0
      };
    }
  }

  private extendCluster(selectedCluster: number): void {
    this.g.selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => {
        if (d.properties.leiden === selectedCluster) {
          // Scale the hexagon coordinates outward
          return this.getScaledPath(d, 1.1); // 10% larger
        }
        // Return original path for non-selected hexagons
        const projection = d3.geoIdentity().fitSize([1200, 1000], {
          type: 'FeatureCollection',
          features: this.features,
        });
        const pathGenerator = d3.geoPath<CellFeature>().projection(projection);
        return pathGenerator(d) || '';
      })
      .style('stroke-width', (d: CellFeature) => {
        return d.properties.leiden === selectedCluster ? '3px' : '1px';
      })
      .style('stroke', (d: CellFeature) => {
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
    let centerX = 0, centerY = 0;
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
      return [
        centerX + dx * scaleFactor,
        centerY + dy * scaleFactor
      ];
    });

    // Create scaled geometry
    const scaledGeometry: CellGeometry = {
      type: 'Polygon',
      coordinates: [scaledCoords]
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

    this.g.selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => pathGenerator(d) || '')
      .style('stroke-width', '1px')
      .style('stroke', 'transparent')
      .style('opacity', 0.8);

    // Reinitialize the mouseleave event
    this.g.selectAll<SVGPathElement, CellFeature>('path')
      .on('mouseleave', (event, d) => this.mouseLeave(event, d));
  }

  private renderNhoodHeatmap(): void {
    if (!this.selectedCell?.properties.leiden_nhood_enrichment) return;

    const enrichment = this.selectedCell.properties.leiden_nhood_enrichment;
    const leiden = this.selectedCell.properties.leiden;
    const n = enrichment.length;
    const clusterLabels = Array.from({ length: n }, (_, i) => `Cluster ${i}`);

    const data: Partial<Plotly.PlotData>[] = [{
      z: [enrichment],
      x: clusterLabels,
      y: [leiden.toString()],
      type: 'heatmap',
      colorscale: 'Viridis'
    }];

    const layout = {
      margin: { t: 30, l: 60, r: 10, b: 40 },
      width: 300,
      height: 170,
      xaxis: {
        title: { text: 'Cluster' },
        automargin: true,
        tickfont: { size: 10 }
      },
      yaxis: {
        title: { text: '' },
        automargin: true,
        showticklabels: false,
        tickfont: { size: 10 }
      }
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
      console.error('Co-occurrence matrix is not an array:', coOccurrenceMatrix);
      this.coOccurrenceData = [];
      return;
    }

    this.coOccurrenceData = [];
    console.log(coOccurrenceMatrix[1][this.selectedInterval])
    try {

      for (let j = 0; j < this.clusterCount; j++) {
        if (Array.isArray(coOccurrenceMatrix[j]) &&
          Array.isArray(coOccurrenceMatrix[j]) &&
          typeof coOccurrenceMatrix[j][this.selectedInterval] === 'number') {
          this.coOccurrenceData.push(coOccurrenceMatrix[j][this.selectedInterval]);
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

    console.log('Co-occurrence data for cluster', this.selectedCluster, 'at interval', this.selectedInterval, ':', this.coOccurrenceData);
  }

  private calculateCoOccurrenceThreshold(): void {
    const allValues = this.coOccurrenceData.flat().filter(val => val > 0);
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
    const blue = Math.floor(255 - (intensity * 200));
    const green = Math.floor(255 - (intensity * 150));
    return `rgb(${blue}, ${green}, 255)`;
  }

  public getIntervalStats(): { min: number, max: number, avg: number } {
    const allValues = this.coOccurrenceData.flat().filter(val => val > 0);
    if (allValues.length === 0) return { min: 0, max: 0, avg: 0 };

    const min = Math.min(...allValues);
    const max = Math.max(...allValues);
    const avg = allValues.reduce((sum, val) => sum + val, 0) / allValues.length;

    return { min, max, avg: Math.round(avg * 100) / 100 };
  }



  private renderLegend(): void {
    // Remove any existing legend
    this.svg.selectAll('.svg-legend').remove();

    if (this.leidenCentralityProps.includes(this.colorByProperty)) {
      // Continuous legend for centrality properties
      const legendX = 20;
      const legendY = 20;
      const width = 250;
      const height = 30;


      const values = this.features.map(f => f.properties.leiden_centrality[this.colorByProperty]);
      const min = Math.min(...values);
      const max = Math.max(...values);

      // Create gradient for continuous legend
      const defs = this.svg.select('defs').empty() ? this.svg.append('defs') : this.svg.select('defs');

      // Remove existing gradient
      defs.select('#svg-legend-gradient').remove();

      const gradient = defs.append('linearGradient')
        .attr('id', 'svg-legend-gradient')
        .attr('x1', '0%')
        .attr('x2', '100%')
        .attr('y1', '0%')
        .attr('y2', '0%');

      // Create gradient stops based on the color scale
      const numStops = 10;
      for (let i = 0; i <= numStops; i++) {
        const t = i / numStops;
        const value = min + t * (max - min);
        gradient.append('stop')
          .attr('offset', `${t * 100}%`)
          .attr('stop-color', this.continuousColorScale(value));
      }

      const legendG = this.svg.append('g')
        .attr('class', 'svg-legend')
        .attr('transform', `translate(${legendX},${legendY})`);

      // Add background for better visibility
      legendG.append('rect')
        .attr('x', -10)
        .attr('y', -25)
        .attr('width', width + 30)
        .attr('height', height + 60)
        .style('fill', 'rgba(255, 255, 255, 0.9)')
        .style('stroke', '#ccc')
        .style('stroke-width', 1)
        .attr('rx', 5);

      // Add the gradient rectangle
      legendG.append('rect')
        .attr('width', width)
        .attr('height', height)
        .style('fill', 'url(#svg-legend-gradient)')
        .style('stroke', '#ccc')
        .style('stroke-width', 1)
        .attr('rx', 3);


      legendG.append('text')
        .attr('x', 0)
        .attr('y', height + 16)
        .attr('text-anchor', 'start')
        .style('font-size', '20px')
        .style('fill', '#333')
        .text(min !== undefined ? min.toFixed(2) : '');


      legendG.append('text')
        .attr('x', width)
        .attr('y', height + 16)
        .attr('text-anchor', 'end')
        .style('font-size', '20px')
        .style('fill', '#333')
        .text(max !== undefined ? max.toFixed(2) : '');


      legendG.append('text')
        .attr('x', width / 2)
        .attr('y', -10)
        .attr('text-anchor', 'middle')
        .style('font-size', '20px')
        .style('font-weight', 'bold')
        .style('fill', '#333')
        .text(this.colorByProperty.replace(/_/g, ' '));

    } else {
      const cellTypes = this.colorScale.domain().sort();
      const legendX = 20;
      const legendY = 6;
      const itemHeight = 30;
      const itemWidth = 200;

      const legendG = this.svg.append('g')
        .attr('class', 'svg-legend')
        .attr('transform', `translate(${legendX},${legendY})`);

      // Add background for categorical legend
      const backgroundHeight = cellTypes.length * itemHeight + 20;
      legendG.append('rect')
        .attr('x', -10)
        .attr('y', -10)
        .attr('width', itemWidth + 20)
        .attr('height', backgroundHeight)
        .style('fill', 'rgba(255, 255, 255, 0.9)')
        .style('stroke', '#ccc')
        .style('stroke-width', 1)
        .attr('rx', 5);

      cellTypes.forEach((cellType, i) => {
        const legendItem = legendG.append('g')
          .attr('transform', `translate(0, ${i * itemHeight})`);

        // Color rectangle
        legendItem.append('rect')
          .attr('width', 30)
          .attr('height', 20)
          .style('fill', this.colorScale(cellType))
          .style('stroke', '#333')
          .style('stroke-width', 0.5)
          .attr('rx', 2);

        // Label
        legendItem.append('text')
          .attr('x', 40)
          .attr('y', 12)
          .style('font-size', '16px')
          .style('fill', '#333')
          .text(cellType);
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

