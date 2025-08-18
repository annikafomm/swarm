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

  public genie3Network: regGraphConnection[] = [];
  public spongeNetwork: regGraphConnection[] = [];
  public geneSetsGenie3: { [regulator: string]: string[] } = {};
public geneSetsSponge: { [regulator: string]: string[] } = {};

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

    // Read sponge network
    d3.csv('assets/sponge_gene_sets_GSM6592049_M2.json', d3.autoType).then((rows) => {
      this.spongeNetwork = rows.map(row => ({
        source: String((row as any)['source'] ?? ''),
        target: String((row as any)['target'] ?? ''),
        weight: Number((row as any)['weight'] ?? 0)
      }));
      console.log('Sponge network loaded:', this.spongeNetwork);
    }).catch((error) => {
      console.error('Error loading sponge network:', error);
    });

    // Read genie3 gene sets
    d3.json<{ [regulator: string]: string[] }>('assets/genie3_gene_sets_GSM6592049_M2.json')
  .then((data) => {
    this.geneSetsGenie3 = data || {};
    console.log('Genie3 gene sets loaded:', Object.keys(this.geneSetsGenie3).length, 'regulators');
  })
  .catch((error) => {
    console.error('Error loading genie3 gene sets:', error);
  });

    // Read sponge gene sets
    d3.json<{ [regulator: string]: string[] }>('assets/sponge_gene_sets_GSM6592049_M2.json')
  .then((data) => {
    this.geneSetsSponge = data || {};
    console.log('Sponge gene sets loaded:', Object.keys(this.geneSetsSponge).length, 'regulators');
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

  public updateAucellGraph(): void {
    // Clear previous graph
    d3.select('#aucell_graph_genie3').selectAll('*').remove();

    if (!this.selectedGeneSetGenie3 || !this.genie3Network) return;

    const regulator = this.selectedGeneSetGenie3;
    const targets = (this.geneSetsGenie3)[regulator] || [];
    // Filter Targets > weight threshold
    const filteredTargets = targets.filter(target => {
      return this.genie3Network.some(connection =>
        connection.source === regulator &&
        connection.target === target &&
        connection.weight > 0.03 // Adjust threshold as needed
      );
    });

    const weightThreshold = 0.03;

    console.log('Regulator:', regulator);
    console.log('Targets:', targets);

    const nodes: { id: string, x?: number, y?: number, group: number }[] = [];
    const edges: { source: string, target: string, weight: number }[] = [];

    // Add regulator node
    nodes.push({ id: regulator, group: 0 });

    // Add target nodes
    filteredTargets.forEach((target: string) => {
      nodes.push({ id: target, group: 1 });
    });

    // Get edges from regulator to targets
    this.genie3Network.forEach(connection => {
      if (connection.source === regulator && filteredTargets.includes(connection.target)) {
        edges.push({
          source: connection.source,
          target: connection.target,
          weight: connection.weight
        });
      }
    });

    // Get neighbors of targets AND regulator with weight > threshold
    const neighborSet = new Set<string>();
    const allMainNodes = [regulator, ...filteredTargets]; // Include regulator in neighbor search

    // Find neighbors for all main nodes (regulator + targets)
    allMainNodes.forEach((mainNode: string) => {
      this.genie3Network.forEach(connection => {

        if (connection.weight > weightThreshold) {

          // If mainNode is the source, add target as neighbor
          if (connection.source === mainNode && !allMainNodes.includes(connection.target)) {
            neighborSet.add(connection.target);
          }
          // If mainNode is the target, add source as neighbor
          else if (connection.target === mainNode && !allMainNodes.includes(connection.source)) {
            neighborSet.add(connection.source);
          }
        }
      });
    });

    console.log('Found neighbors:', Array.from(neighborSet));
    console.log('Neighbor count:', neighborSet.size);

    // Add neighbor nodes
    if (nodes.length < 30) {
      Array.from(neighborSet).forEach(neighbor => {
        nodes.push({ id: neighbor, group: 2 });
      });
    }

    // Add ALL connections between any nodes in our network with sufficient weight
    const allNodeIds = new Set(nodes.map(n => n.id));

    this.genie3Network.forEach(connection => {
      if (connection.weight > weightThreshold) {
        const sourceInNetwork = allNodeIds.has(connection.source);
        const targetInNetwork = allNodeIds.has(connection.target);

        // Add edge if both nodes are in our network
        if (sourceInNetwork && targetInNetwork) {
          // Avoid duplicate edges
          const edgeExists = edges.some(e =>
            (e.source === connection.source && e.target === connection.target) ||
            (e.source === connection.target && e.target === connection.source)
          );

          if (!edgeExists) {
            edges.push({
              source: connection.source,
              target: connection.target,
              weight: connection.weight
            });
          }
        }
      }
    });

    console.log('Final nodes:', nodes.length);
    console.log('Final edges:', edges.length);
    console.log('Edges:', edges);



    console.log('After enhancement - Nodes:', nodes.length, 'Edges:', edges.length);


    if (nodes.length === 0) {
      console.warn('No nodes to display');
      return;
    }


    // Create the graph
    const graph = {
      nodes: nodes.filter(node => node.id && node.id.length > 0),
      edges: edges.filter(edge =>
        nodes.some(node => node.id === edge.source) &&
        nodes.some(node => node.id === edge.target)
      )
    };

    // Create the graph visualization
    const width = 500;
    const height = 300;
    const margin = 20;

    const svg = d3.select('#aucell_graph_genie3')
      .append('svg')
      .attr('width', width)
      .attr('height', height)
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
        switch(d.group) {
          case 0: return 15; // regulator
          case 1: return 12; // targets
          case 2: return 8;  // neighbors
          case 3: return 6;
          default: return 10;
        }
      })
      .attr('fill', (d: any) => {
        switch(d.group) {
          case 0: return '#e41a1c'; // regulator - red
          case 1: return '#377eb8'; // targets - blue
          case 2: return '#4daf4a'; // neighbors - green
          case 3: return '#ff7f00'; // high-weight - orange
          default: return '#999';
        }
      });

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

    // Add weight labels on edges
    const edgeLabels = svg.append('g')
      .selectAll('text')
      .data(graph.edges)
      .enter()
      .append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '0.35em')
      .style('font-size', '8px')
      .style('fill', '#666')
      .style('font-weight', 'bold')
      .style('background-color', 'white')
      .style('padding', '1px')
      .text((d: any) => d.weight.toFixed(2));

    // Initialize simulation with stronger forces
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

      // Position edge labels at the midpoint of each edge
      edgeLabels
        .attr('x', (d: any) => (d.source.x + d.target.x) / 2)
        .attr('y', (d: any) => (d.source.y + d.target.y) / 2);
    });

    console.log('Network visualization complete');
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

  public openSidenav(event: MouseEvent, cell: CellFeature): void {
    this.selectedCell = cell;
    // First get the associated gene sets
    [this.selectedCellAssociatedGeneSetsGenie3, this.selectedCellAssociatedGeneSetsSponge] = this.getAssociatedGeneSets(cell);

    // Then set the selected gene set if we have any
    if (this.selectedCellAssociatedGeneSetsGenie3.length > 0) {
      this.selectedGeneSetGenie3 = this.selectedCellAssociatedGeneSetsGenie3[0]; // Select the first one
    } else {
      this.selectedGeneSetGenie3 = null; // Clear if no gene sets found
    }

    if (this.colorByProperty === 'leiden') {
      this.openClusterSidenav(cell.properties.leiden);
      this.extendCluster(cell.properties.leiden);
    }
    else {
      d3.select(event.target as SVGElement)
        .transition()
        .style('stroke', 'black');
    }


    setTimeout(() => this.renderNhoodHeatmap(), 0);

    setTimeout(() => this.updateAucellGraph(), 0);
  }

  public openClusterSidenav(clusterId: number): void {
    this.selectedCluster = clusterId;
    this.clusterCells = this.features.filter(cell => cell.properties.leiden === clusterId);
    this.calculateClusterStats();

    // Initialize co-occurrence table for this cluster
    this.updateCoOccurrenceTable();

    if (this.clusterCells.length > 0) {
      this.selectedCell = this.clusterCells[0];
      setTimeout(() => this.renderNhoodHeatmap(), 100);
      setTimeout(() => this.updateAucellGraph(), 100);
    }

  }

  public onGeneSetChange(): void {
  console.log('Gene set changed to:', this.selectedGeneSetGenie3);

  // Clear previous graph and regenerate
  if (this.selectedGeneSetGenie3) {
    setTimeout(() => this.updateAucellGraph(), 0);
  }
}

  public selectCellFromCluster(cell: CellFeature): void {
    this.selectedCell = cell;
    setTimeout(() => this.renderNhoodHeatmap(), 0);
    setTimeout(() => this.updateAucellGraph(), 0);
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

interface regGraphConnection {
  source: string;
  target: string;
  weight: number;
}
