import { Component, ElementRef, HostListener } from '@angular/core';
import { CommonModule } from '@angular/common';
import { GeoDataService } from '../geo-data.service';
import { HttpClient } from '@angular/common/http';
import { SessionService } from '../session.service';

@Component({
  selector: 'app-download-menu',
  imports: [CommonModule],
  templateUrl: './download-menu.component.html',
  styleUrl: './download-menu.component.scss',
})
export class DownloadMenuComponent {
  menuOpen = false;

  constructor(
    private eRef: ElementRef,
    private geoDataService: GeoDataService,
    private http: HttpClient,
    private sessionService: SessionService,
  ) {}

  toggleMenu() {
    this.menuOpen = !this.menuOpen;
  }

  @HostListener('document:click', ['$event'])
  onClickOutside(event: Event) {
    if (this.menuOpen && !this.eRef.nativeElement.contains(event.target)) {
      this.menuOpen = false;
    }
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

    const data = this.geoDataService.getData();
    if (!data) {
      console.warn('No GeoJSON data available for download');
      return;
    }

    const jsonString = JSON.stringify(data, null, 4);
    const jsonBlob = new Blob([jsonString], { type: 'application/json' });
    const jsonUrl = URL.createObjectURL(jsonBlob);

    const link = document.createElement('a');
    link.href = jsonUrl;
    link.download = 'map.geojson';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);

    URL.revokeObjectURL(jsonUrl);
  }

  downloadAnndata() {
    this.menuOpen = false;

    this.sessionService
      .callWithSession(() =>
        this.http.get(`${this.sessionService.apiUrl}/download_adata`, {
          withCredentials: true,
          responseType: 'blob',
        }),
      )
      .subscribe((blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = 'adata.h5ad';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        window.URL.revokeObjectURL(url);
      });
  }
}
