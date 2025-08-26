import { Component, ViewChild, ElementRef } from '@angular/core';
import { CommonModule } from '@angular/common';

@Component({
  selector: 'app-download-menu',
  imports: [CommonModule],
  templateUrl: './download-menu.component.html',
  styleUrl: './download-menu.component.scss',
})
export class DownloadMenuComponent {
  //@ViewChild('#hexbin') hexbinRef!: ElementRef<SVGSVGElement>;
  @ViewChild('hexbin', { static: false })
  hexbinRef!: ElementRef<HTMLDivElement>;

  //@ViewChild('legend') legendRef!: ElementRef<SVGGElement>;

  menuOpen = false;

  toggleMenu() {
    this.menuOpen = !this.menuOpen;
  }

  doAction(action: string) {
    this.menuOpen = false; // close menu after action
    switch (action) {
      case 'one':
        console.log('Action One executed!');
        break;
      case 'two':
        alert('Action Two executed!');
        break;
      case 'three':
        console.log('Action Three executed! Performing calculation:', 2 + 2);
        break;
      case 'four':
        window.open('https://angular.dev', '_blank');
        break;
    }
  }

  downloadSvg() {
    this.menuOpen = false;
    const svg = document.querySelector<HTMLDivElement>('#hexbin > svg');
    if (!svg) {
      console.warn('No <svg> found inside #hexbin');
      return;
    }

    const clonedSvg = svg.cloneNode(true) as SVGSVGElement;

    const svgData = new XMLSerializer().serializeToString(clonedSvg);
    const svgBlob = new Blob([svgData], {
      type: 'image/svg+xml;charset=utf-8',
    });
    const svgUrl = URL.createObjectURL(svgBlob);

    const link = document.createElement('a');
    link.href = svgUrl;
    link.download = 'map.svg';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);

    URL.revokeObjectURL(svgUrl);
  }

  downloadPng() {
    this.menuOpen = false;
    // Grab the D3-generated SVG inside #hexbin
    const svgElement = document.querySelector<SVGSVGElement>('#hexbin > svg');
    if (!svgElement) {
      console.warn('No SVG found inside #hexbin');
      return;
    }

    const width = 1500;
    const height = 1200;

    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    // Clone the SVG to avoid mutating the original
    const clonedSvg = svgElement.cloneNode(true) as SVGSVGElement;

    // Serialize and convert to data URL
    const svgData = new XMLSerializer().serializeToString(clonedSvg);
    const svgDataUrl =
      'data:image/svg+xml;base64,' +
      btoa(unescape(encodeURIComponent(svgData)));

    const img = new Image();
    img.onload = () => {
      ctx.clearRect(0, 0, width, height);
      ctx.drawImage(img, 0, 0, width, height);

      // Trigger download
      const pngDataUrl = canvas.toDataURL('image/png');
      const link = document.createElement('a');
      link.href = pngDataUrl;
      link.download = 'map.png';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    };

    img.src = svgDataUrl;
  }

  downloadGeojson() {
    this.menuOpen = false;
    alert('TODO');
  }

  downloadAnndata() {
    this.menuOpen = false;
    alert('TODO');
  }
}
