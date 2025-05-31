import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import * as d3 from 'd3';

@Component({
  selector: 'app-hexagon-plot',
  imports: [CommonModule,d3],
  standalone: true,
  templateUrl: './hexagon-plot.component.html',
  styleUrl: './hexagon-plot.component.scss'
})
export class HexagonPlotComponent {
  public hexagonData: hexMap[] = [];
  private hexagonSize: number = 7;
  public selectedHexagon: hexMap | null = null;

  constructor() {
    this.generateRandomCoordinates();
   } // You have to play with these values to center your map

   openSidenav(hex: hexMap) {
    this.selectedHexagon = hex;
  }

  closeSidenav() {
    this.selectedHexagon = null;
  }

  generateRandomCoordinates() {

    this.hexagonData = [];
    const cols = 100; // number of hexagons per row
    const rows = 100;  // number of rows
    const size = this.hexagonSize;
    const w = size * 1.5;
    const h = size * Math.sqrt(3);

    let i = 0;
    for (let row = 0; row < rows; row++) {
      for (let col = 0; col < cols; col++) {
        // Offset every other row
        const x = col * w + (row % 2) * (w / 2) + 30;
        const y = row * h + 30;
        const hexagon = {
          type: 'FeatureCollection',
          geometry: {
            type: 'Polygon',
            coordinates: [
              [
                [x, y],
                [x + size, y],
                [x + size * 1.5, y + size * Math.sqrt(3) / 2],
                [x + size, y + size * Math.sqrt(3)],
                [x, y + size * Math.sqrt(3)],
                [x - size * 0.5, y + size * Math.sqrt(3) / 2],
                [x, y]
              ]
            ]
          },
          properties: {
            id: i,
            label: `Hexagon ${i}`,
            color: d3.interpolateRainbow(i / (rows * cols))
          }
        };
        this.hexagonData.push(hexagon);
        i++;
    }
  }
}

}

export interface hexMap {
  type: string;
  geometry: {
    type: string;
    coordinates: number[][][];
  };
  properties: {
    id: number;
    label: string;
    color: string;
  }
}
