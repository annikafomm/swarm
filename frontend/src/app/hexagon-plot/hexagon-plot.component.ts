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

  private colorScale = d3
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


        const projection = d3.geoIdentity().fitSize([width,height], {
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
          .attr('d', (d) => pathGenerator(d)) // Fixed: Now properly typed
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
    if (this.leidenCentralityProps.includes(this.colorByProperty)) {
      // Get all values for the selected centrality property
      const values = this.features.map(f => f.properties.leiden_centrality[this.colorByProperty]);
      const min = Math.min(...values);
      const max = Math.max(...values);

      // Set domain for continuous scale
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

     d3.select(event.target as SVGElement)
      .transition()
      .style('stroke', 'black');

    this.selectedCell = cell;
    setTimeout(() => this.renderNhoodHeatmap(), 0); // Wait for DOM update

    setTimeout(() => this.updateAucellGraph(), 0);
  }

  private renderNhoodHeatmap(): void {
    if (!this.selectedCell?.properties.leiden_nhood_enrichment) return;

    const enrichment = this.selectedCell.properties.leiden_nhood_enrichment;
    const leiden = this.selectedCell.properties.leiden;
    const n = enrichment.length;
    const clusterLabels = Array.from({ length: n }, (_, i) => `Cluster ${i}`);

    // 1D heatmap: just a single row
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

    const container = document.getElementById('nhood-heatmap');
    if (!container) return;
    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, { displayModeBar: false });
  }



  public closeSidenav(cell:CellFeature): void {
    this.selectedCell = null;
    this.updateHexColors();
  }



  private renderLegend(): void {
    // Remove previous legend group
    this.svg.selectAll('.svg-legend').remove();

    const svgWidth = +this.svg.attr('width');
    const svgHeight = +this.svg.attr('height');

    if (this.leidenCentralityProps.includes(this.colorByProperty)) {
      // Continuous legend
      const width = 160, height = 14;
      const legendX = svgWidth - width;
      const legendY = svgHeight - height - 38;

      const legendG = this.svg.append('g')
        .attr('class', 'svg-legend')
        .attr('transform', `translate(${legendX},${legendY})`);

      // Gradient definition
      const defs = this.svg.select('defs').empty()
        ? this.svg.append('defs')
        : this.svg.select('defs');
      defs.select('linearGradient#svg-legend-gradient').remove();
      const gradient = defs.append('linearGradient')
        .attr('id', 'svg-legend-gradient')
        .attr('x1', '0%').attr('y1', '0%')
        .attr('x2', '100%').attr('y2', '0%');

      const [min, max] = this.continuousColorScale.domain();
      for (let i = 0; i <= 100; i++) {
        gradient.append('stop')
          .attr('offset', `${i}%`)
          .attr('stop-color', this.continuousColorScale(min + (i / 100) * (max - min)));
      }

      legendG.append('rect')
        .attr('width', width)
        .attr('height', height)
        .style('fill', 'url(#svg-legend-gradient)')
        .style('stroke', '#ccc')
        .style('stroke-width', 1)
        .attr('rx', 5);

      legendG.append('text')
        .attr('x', 0)
        .attr('y', height + 14)
        .attr('text-anchor', 'start')
        .style('font-size', '11px')
        .text(min !== undefined ? min.toFixed(2) : '');

      legendG.append('text')
        .attr('x', width)
        .attr('y', height + 14)
        .attr('text-anchor', 'end')
        .style('font-size', '11px')
        .text(max !== undefined ? max.toFixed(2) : '');

      legendG.append('text')
        .attr('x', width / 2)
        .attr('y', -4)
        .attr('text-anchor', 'middle')
        .style('font-size', '12px')
        .style('font-weight', 'bold')
        .text(this.colorByProperty.replace(/_/g, ' '));
    } else {
      // Categorical legend
      const cellTypes = this.colorScale.domain().sort();
      const legendX = 0;
      const legendY = 60;
      const width = 100;
      const itemHeight = 25;

      const legendG = this.svg.append('g')
        .attr('class', 'svg-legend')
        .attr('transform', `translate(${legendX},${legendY})`);

      cellTypes.forEach((d, i) => {
        legendG.append('rect')
          .attr('x', 0)
          .attr('y', i * itemHeight + 4)
          .attr('width', 24)
          .attr('height', 20)
          .attr('fill', this.colorScale(d));

        legendG.append('text')
          .attr('x', 32)
          .attr('y', i * itemHeight + 18)
          .text(d)
          .style('font-size', '16px')
          .attr('alignment-baseline', 'middle');
      });

      legendG.append('text')
        .attr('x', width / 2)
        .attr('y', -6)
        .attr('text-anchor', 'middle')
        .style('font-size', '12px')
        .style('font-weight', 'bold')
        .text(this.colorByProperty.replace(/_/g, ' '));
    }
  }


}

interface CellGeometry {
  type: 'Polygon';
  coordinates: number[][][];
}

interface CellProperties {
  barcode: string;
  centroid: [];
  cell_type: string;
  leiden_nhood_enrichment: number[];
  leiden: number;
  color: string;
  aucell: { [key: string]: number };
  leiden_centrality: { [key: string]: number };
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
