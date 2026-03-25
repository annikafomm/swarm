import { Component, OnInit, OnDestroy, HostListener, HostBinding } from '@angular/core';
import { RouterOutlet, RouterLink } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
import { DownloadMenuComponent } from './download-menu/download-menu.component';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';
import { PathsService } from './paths.service';
import { DatasetService } from './datasets.service';
import { MatDialog } from '@angular/material/dialog';
import { UnregisteredDatasetsDialogComponent } from './unregistered-datasets-dialog/unregistered-datasets-dialog.component';
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
    private datasetService: DatasetService,
    private dialog: MatDialog,
  ) { }

  ngOnInit() {
    this.sessionService.initSession();

    // Show unregistered datasets dialog after session is initialized
    // Track by session ID - show dialog on each NEW session (page refresh creates new session)
    setTimeout(() => {
      const currentSessionId = this.sessionService.sessionId;
      const lastDialogSessionId = sessionStorage.getItem('lastUnregisteredDialogSessionId');

      console.log(`[APP] Current session: ${currentSessionId}, Last dialog session: ${lastDialogSessionId}`);

      // Show dialog if this is a new session (different from last time we showed it)
      if (currentSessionId && currentSessionId !== lastDialogSessionId) {
        console.log('[APP] New session detected, showing unregistered datasets dialog');
        this.showUnregisteredDatasetsDialog();
        sessionStorage.setItem('lastUnregisteredDialogSessionId', currentSessionId);
      } else {
        console.log('[APP] Same session, skipping unregistered datasets dialog');
      }
    }, 1000);  // Increased timeout to ensure session is ready

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

  private showUnregisteredDatasetsDialog(): void {
    console.log('[DEBUG] Checking for unregistered datasets...');
    this.sessionService.callWithSession(() =>
      this.datasetService.loadUnregisteredDatasets()
    )
      .pipe(takeUntil(this.destroy$))
      .subscribe({
        next: (response) => {
          console.log('[DEBUG] Got response:', response);
          const datasets = response.datasets || [];
          console.log(`[DEBUG] Found ${datasets.length} unregistered datasets`);
          // Always show dialog on new session, regardless of dataset count
          console.log('[DEBUG] Opening unregistered datasets dialog on new session');
          const dialogRef = this.dialog.open(UnregisteredDatasetsDialogComponent, {
            width: '1000px',
            maxHeight: '80vh',
            disableClose: false,
            autoFocus: 'first-button',
          });

          // Reload datasets when dialog closes
          dialogRef.afterClosed()
            .pipe(takeUntil(this.destroy$))
            .subscribe((result) => {
              if (result?.datasetsChanged) {
                // Reload available datasets to reflect any registrations/deletions
                this.datasetService.loadAvailableDatasets();
                console.log('✓ Datasets list refreshed after dataset management');
              }
            });
        },
        error: (error) => {
          console.error('[DEBUG] Error loading unregistered datasets:', error);
          console.warn('Could not check for unregistered datasets:', error);
          // Silently fail - not critical to functionality
        }
      });
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
