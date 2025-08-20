import { Component, OnInit, ɵisComponentDefPendingResolution } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ActivatedRoute } from '@angular/router';
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
  public Math = Math;
  public selectedCell: CellFeature | null = null;
  private svg!: d3.Selection<SVGSVGElement, unknown, HTMLElement, any>;
  private g!: d3.Selection<SVGGElement, unknown, HTMLElement, any>;
  public datasetTitle = 'GSM6592049_M2';
  public genie3Network: genie3Connection[] = [];

  public colorableProperties = ['cell_type', 'leiden', 'degree_centrality', 'average_clustering', 'closeness_centrality'];
  public leidenCentralityProps = [
    'degree_centrality',
    'average_clustering',
    'closeness_centrality'
  ];

  public colorByProperty = 'cell_type';
  public selectedGeneSet: string | null = null;

  public selectedCluster: number | null = null;
  public clusterCells: CellFeature[] = [];
  public clusterCellTypes: { type: string, count: number, percentage: string }[] = [];
  public clusterCentralityAvg: { degree_centrality: number, average_clustering: number, closeness_centrality: number } = {
    degree_centrality: 0,
    average_clustering: 0,
    closeness_centrality: 0
  };

  

  public selectedInterval: number = 0;
  public coOccurrenceData: number[] = [];
  public coOccurrenceColumns: string[] = [];
  public coOccurrenceThreshold: number = 0.5;
  public maxInterval: number = 49;
  public clusterCount: number = 10;

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
  private features: CellFeature[] = [];

  constructor() {

  }
  ngOnInit(): void {
    this.createHexagonPlot();
    this.loadAndRenderData();

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

  private loadAndRenderData(): void {
    d3.json<GeoJsonData>('assets/hexagons_GSM6592049_M2.geojson')
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
    d3.select('#aucell_graph').selectAll('*').remove();
    if (!this.selectedGeneSet || !this.genie3Network) return;

    const regulator = this.selectedGeneSet;
    const edges = this.genie3Network.filter(conn => conn.source === regulator && conn.weight > 0.02);

    if (edges.length === 0) return;

    const width = 260, height = 180, centerX = width / 2, centerY = height / 2, r = 60;
    const n = edges.length;

    // Central node (regulator)
    const nodes = [{ id: regulator, x: centerX, y: centerY, group: 1 }];
    // Place targets in a circle
    edges.forEach((e, i) => {
      const angle = (2 * Math.PI * i) / n;
      nodes.push({
        id: e.target,
        x: centerX + r * Math.cos(angle),
        y: centerY + r * Math.sin(angle),
        group: 2
      });
    });

    const svg = d3.select('#aucell_graph')
      .append('svg')
      .attr('width', width)
      .attr('height', height);

    // Draw edges
    svg.selectAll('line')
      .data(edges)
      .enter()
      .append('line')
      .attr('x1', centerX)
      .attr('y1', centerY)
      .attr('x2', (_d, i) => nodes[i + 1].x)
      .attr('y2', (_d, i) => nodes[i + 1].y)
      .attr('stroke', '#888')
      .attr('stroke-width', d => Math.max(1, d.weight * 2));

    // Draw nodes
    svg.selectAll('circle')
      .data(nodes)
      .enter()
      .append('circle')
      .attr('cx', d => d.x)
      .attr('cy', d => d.y)
      .attr('r', d => d.group === 1 ? 18 : 12)
      .attr('fill', d => d.group === 1 ? '#FF7373' : '#66cdaa')
      .attr('stroke', '#333')
      .attr('stroke-width', 1.5);

    // Draw labels
    svg.selectAll('text')
      .data(nodes)
      .enter()
      .append('text')
      .attr('x', d => d.x)
      .attr('y', d => d.y + (d.group === 1 ? 5 : 4))
      .attr('text-anchor', 'middle')
      .style('font-size', d => d.group === 1 ? '13px' : '11px')
      .text(d => d.id);
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
  aucell: { [key: string]: number };
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

interface genie3Connection {
  source: string;
  target: string;
  weight: number;
}
