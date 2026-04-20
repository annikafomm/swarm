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
  public unregisteredDialogOpen = false; // Track if unregistered datasets dialog is open

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
        if (paths.adataPath) promises.push(loadPath('read_adata', paths.adataPath, 'adata'));
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
    this.sessionService.callWithSession(() =>
      this.datasetService.loadUnregisteredDatasets()
    )
      .pipe(takeUntil(this.destroy$))
      .subscribe({
        next: (response) => {
          const datasets = response.datasets || [];
          console.log(`[DEBUG] Found ${datasets.length} unregistered datasets`);
          this.unregisteredDialogOpen = true;
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
              this.unregisteredDialogOpen = false;
              if (result?.datasetsChanged) {
                // Reload available datasets to reflect any registrations/deletions
                this.datasetService.loadAvailableDatasets();
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
      text: 'Welcome to SWARM! I am your helper bee, here to guide you through the main features of this tool. Let\'s take a quick tour to get you started.',
      buttons: [
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
      ]
    });

    tour.addStep({
      id: 'main-view',
      text:
        'This is the <b>main view</b> visualizing the spatial coordinates of your data as <b>hexagons</b>. It allows you to explore your data interactively.<ul><li>All <b>obs columns</b> and computed features are available for coloring the hexagons (if it is cell/spot level information)</li><li>The <b>zoom and pan</b> functionality allows you to navigate through the plot</li><li>Clicking on a hexagon allows you to display detailed information about the spot/cell in the <b>Cell Information Tab</b></li><li>The <b>legend</b> on the left helps you interpret the colors on the map (continuous values are normalized between 0 and 1, categorical values are assigned distinct colors)</li></ul>',
      attachTo: { element: '#hexbin', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
      ]
    });

    tour.addStep({
      id: 'upload-data-btn',
      text:
        'To <b>upload your data</b>, click on the <b>"Upload data"</b> button in the top navigation bar.<br>To get a detailed tutorial for the upload process, click on the <b>"Upload"</b> button.',
      attachTo: { element: '#upload-data-btn', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        {
          text: 'Upload', action: () => {
            tour.complete();
            this.uploadTutorial();
            // Handle tutorial button click
          }
        },
        { text: 'Close', action: tour.complete }
      ]
    });


    tour.addStep({
      id: 'selectViews',
      text:
        'To explore the different parts of your analysis, as defined in the upload forms, there are several different <b>modes</b> available.<br>Use the <b>tabs</b> to switch between modes such as <b>gene expression</b>, <b>regulatory scores</b>, <b>LIANA+ scores</b>, <b>DGEA</b>, and <b>Multiome</b> specific modes.<br>Each mode provides unique insights into the regulatory landscape of your data.',
      attachTo: { element: '.metadata-tables', on: 'left' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
      ]
    });

    tour.addStep({
      id: 'dataset-selection',
      text:
        'If you have uploaded multiple datasets, you can select which one to visualize using the <b>dataset selection dropdown</b>.<br><b>Previous uploads</b> can be recovered in the pop-up dialog that appears on initialization of the tool.'
      ,
      attachTo: { element: '#dataset-dropdown-main', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
      ]
    });

    tour.addStep({
      id: 'compare-view',
      text:
        'If there is something within or between datasets you want to directly compare, this is possible by using the comparison feature.' +
        'You can toggle the comparison view using this button.',
      attachTo: { element: '#compare-button', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
      ]
    });

    tour.addStep({
      id: 'info-page',
      text:
        'For more information about SWARM, including documentation and support resources, click on the "Info" button in the top navigation bar.',
      attachTo: { element: '.top-bar .info-btn', on: 'bottom' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        { text: 'Done', action: tour.complete }
      ]
    });

    tour.addStep({
      id: 'download-menu',
      text:
        'The <b>download menu</b> allows you to explore export options for your data.<ul><li>Export files of the current visualization (<b>AnnData, GeoJSON, SVG</b>, etc.)</li><li>Select your dataset and go to the <b>Info-Page</b> to access other files from your analysis</li></ul>',
      attachTo: { element: 'app-download-menu', on: 'left' },
      buttons: [
        { text: 'Back', action: tour.back },
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
      ]
    });



    tour.start();
  }

  public uploadTutorial(): void {
    // Open forms
    this.openForm();

    const tour = new Shepherd.Tour({
      useModalOverlay: true,
      defaultStepOptions: {
        classes: 'shepherd-theme-custom',
        scrollTo: true
      }
    });

    tour.addStep({
      id: 'upload-tutorial-intro',
      text: 'This tutorial will guide you through the process of uploading your data to SWARM. We will cover the different types of data you can upload and how to fill out the upload form correctly.',
      buttons: [
        { text: 'Next', action: tour.next },
        { text: 'Close', action: tour.complete }
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
