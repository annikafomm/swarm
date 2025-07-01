import { Component, OnInit } from '@angular/core';
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
  private g!: d3.Selection<SVGGElement, unknown, HTMLElement,any>;
  public datasetTitle = 'GSM6592049_M2';


  public colorableProperties = ['cell_type', 'leiden'];
  public colorByProperty = 'cell_type';

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

  private features: CellFeature[] = [];

  constructor() {

  }
  ngOnInit(): void {
    this.createHexagonPlot();
    this.loadAndRenderData();

  }
  private createHexagonPlot(): void {
    this.svg = d3
      .select('#hexbin')
      .append('svg')
      .attr('width', 700)
      .attr('height', 500)
      .attr('viewBox', [0, 0, 1000, 700] as [number, number, number, number])
      .style('background-color', 'white')
      .style('overflow', 'hidden')
      .call(
        d3
          .zoom<SVGSVGElement, unknown>()
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


        const projection = d3.geoIdentity().fitSize([1000, 800], {
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
            const total = d.properties.cell_type|| 0;
            return this.colorScale(total.toString());
          })
          .style('stroke', 'transparent')
          .style('opacity', 0.8)
          .on('mouseover', (event, d) => this.mouseOver(event, d))
          .on('mouseleave', (event, d) => this.mouseLeave(event, d))
          .on('click', (event,d) => this.openSidenav(event,d))

          this.colorScale.domain([...new Set(data.features.map((f: CellFeature) => f.properties.cell_type))]);
          this.renderLegend();
      })
      .catch((error) => {
        console.error('Error loading or rendering data:', error);
      });



  }

   public updateHexColors(): void {
    // Update color scale domain for the new property
    this.colorScale.domain([...new Set(this.features.map(f => String(f.properties[this.colorByProperty])))]);
    // Update hexagon colors
    this.g.selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('fill', (d: CellFeature) => this.colorScale(String(d.properties[this.colorByProperty])));

    this.renderLegend();
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
    setTimeout(() => this.renderNhoodHeatmap(), 0); // Wait for DOM update
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



  public closeSidenav(): void {
    this.selectedCell = null;
  }

  private renderLegend(): void {
    // Clear existing legend
    this.svg.selectAll('svg').remove();
    // Get unique cell types from your color scale's domain
    const cellTypes = this.colorScale.domain().sort();


    const legendItem = this.svg
      .append('svg')
      .attr('width', 400)
      .attr('height', cellTypes.length * 28)
      .selectAll('g')
      .data(cellTypes)
      .enter()
      .append('g')
      .attr('transform', (_d, i) => `translate(0,${i * 28})`);

    legendItem
      .append('rect')
      .attr('x', 0)
      .attr('y', 4)
      .attr('width', 24)
      .attr('height', 20)
      .attr('fill', d => this.colorScale(d));

    legendItem
      .append('text')
      .attr('x', 32)
      .attr('y', 18)
      .text(d => d)
      .style('font-size', '16px')
      .attr('alignment-baseline', 'middle');
  }
}

interface CellGeometry {
  type: 'Polygon';
  coordinates: number[][][];
}

interface CellProperties {
  barcode : string;
  centroid: [];
  cell_type: string;
  leiden_nhood_enrichment: number[];
  leiden: number;
  color: string;
  [key: string]: string | number | number[] | [] | undefined;
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
