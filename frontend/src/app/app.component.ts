import { Component, OnInit, OnDestroy, HostListener, HostBinding } from '@angular/core';
import { RouterOutlet, RouterLink } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
import { DownloadMenuComponent } from './download-menu/download-menu.component';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';
import { PathsService } from './paths.service';
import { Subject } from 'rxjs';
import { takeUntil } from 'rxjs/operators';
import Shepherd from 'shepherd.js';

@Component({
  selector: 'app-root',
  standalone: true,
  imports: [
    CommonModule,
    FormPageComponent,
    DownloadMenuComponent,
    RouterOutlet,
    RouterLink,
  ],
  templateUrl: './app.component.html',
  styleUrl: './app.component.scss',
})
export class AppComponent implements OnInit, OnDestroy {
  private destroy$ = new Subject<void>();
  private isLoadingPaths = false; // Prevent concurrent path loads

  constructor(
    private http: HttpClient,
    public sessionService: SessionService,
    private pathsService: PathsService,
  ) { }

  ngOnInit() {
    this.sessionService.initSession();

    this.pathsService.paths$
      .pipe(takeUntil(this.destroy$))
      .subscribe((paths) => {
        if (!paths || this.isLoadingPaths) return;

        this.isLoadingPaths = true;

        const loadPath = (backendUrl: string, path: string, label: string) => {
          return this.sessionService.callWithSession(() =>
            this.http.post(
              `${this.sessionService.apiUrl}/${backendUrl}`,
              { path },
              { withCredentials: true }
            )
          ).toPromise().then(
            res => {
              console.log(`[Backend] Loaded ${label}`, res);
              return res;
            },
            err => {
              console.error(`[Backend] Failed to load ${label} ${path}`, err);
              return null;
            }
          );
        };

        // Load paths sequentially to avoid race conditions
        const promises: Promise<any>[] = [];
        if (paths.adataMainPath) promises.push(loadPath('read_adata', paths.adataMainPath, 'adata'));
        if (paths.genieFiltPath) promises.push(loadPath('read_network_genie', paths.genieFiltPath, 'network_genie'));
        if (paths.spongeFiltPath) promises.push(loadPath('read_network_sponge', paths.spongeFiltPath, 'network_sponge'));

        Promise.all(promises).finally(() => {
          this.isLoadingPaths = false;
        });
      });
  }

  ngOnDestroy() {
    this.destroy$.next();
    this.destroy$.complete();
  }

  title = 'frontend';
  public formsOpen = false;

  @HostBinding('class.modal-open')
  get modalOpen(): boolean {
    return this.formsOpen;
  }

  public closeSidenavForms(): void {
    this.formsOpen = false;
  }

  public openForm(): void {
    this.formsOpen = true;
  }

  public basicTutorial(): void {
    const tour = new Shepherd.Tour({
      useModalOverlay: true,
      defaultStepOptions: {
        classes: 'shepherd-theme-custom',
        scrollTo: true
      }
    });

    tour.addStep({
      id: 'hexagon-plot-intro',
      text: 'Welcome to SWARM! In this tutorial, we will guide you through the main features of the tool.',
      buttons: [
        { text: 'Next', action: tour.next }
      ]
    });

    tour.addStep({
      id: 'main-view',
      text:
        'This is the main view visualizing the spatial transcriptomics data as a hexagon plot. It allows you to explore your data interactively.' +
        ' Depending on the selected view, different properties of the cells are represented by the colors in the graph as also indicated by the legend in the top left corner.' +
        ' You can click individual hexagons to get more information about the cells or spots they represent.' +
        ' Additionally you can use the zoom and pan functionality to navigate through the plot.',
      attachTo: { element: '#hexbin', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next }
      ]
    });

    tour.addStep({
      id: 'selectViews',
      text:
        'You can select different views to visualize various aspects of your data.' +
        ' Use the tabs to switch between views such as gene expression, cell type, regulatory scores and LIANA+ scores.' +
        ' Each view provides unique insights into the spatial organization and characteristics of your data.',
      attachTo: { element: '.metadata-tables', on: 'left' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next }
      ]
    });

    tour.addStep({
      id: 'download-menu',
      text:
        'The download menu allows you to export the current visualization in various formats such as SVG, PNG, GeoJSON, and AnnData.' +
        ' This enables you to save your results and share them with others or use them for further analysis.',
      attachTo: { element: 'app-download-menu', on: 'left' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next }
      ]
    });

    tour.addStep({
      id: 'upload-data',
      text:
        'To upload new data, click on the "Upload data" button in the top navigation bar.' +
        ' This will open a form where you can select and upload your spatial transcriptomics data files.' +
        ' Follow the instructions in the form to ensure successful data upload and processing.',
      attachTo: { element: '.top-bar button', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next }
      ]
    });

    tour.addStep({
      id: 'info-page',
      text:
        'For more information about SWARM, including documentation and support resources, click on the "Info" button in the top navigation bar.' +
        ' This will take you to the info page where you can find helpful materials to assist you in using the tool effectively.',
      attachTo: { element: '.top-bar .info-btn', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Done', action: tour.complete }
      ]
    });

    tour.start();
  }

  @HostListener('window:beforeunload', ['$event'])
  beforeUnloadHandler(_event: Event) {
    // Use sendBeacon for reliable delivery during page unload
    // It works even if the page is being closed/refreshed
    const sessionId = this.sessionService.sessionId;
    if (sessionId) {
      const url = `${this.sessionService.apiUrl}/delete_session`;
      navigator.sendBeacon(url, JSON.stringify({ session_id: sessionId }));
      console.log('[Frontend] Sent delete_session beacon');
    }
  }
}
