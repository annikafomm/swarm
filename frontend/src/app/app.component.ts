import { Component, OnInit, HostListener } from '@angular/core';
import { RouterOutlet, RouterLink } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
import { DownloadMenuComponent } from './download-menu/download-menu.component';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';
import { PathsService } from './paths.service';

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
export class AppComponent implements OnInit {
  constructor(
    private http: HttpClient,
    public sessionService: SessionService,
    private pathsService: PathsService,
  ) {}

  ngOnInit() {
    this.sessionService.initSession();

    this.pathsService.paths$.subscribe((paths) => {
      if (!paths) return;

      const loadPath = (backendUrl: string, path: string, label: string) => {
        this.sessionService.callWithSession(() =>
          this.http.post(
            `${this.sessionService.apiUrl}/${backendUrl}`,
            { path },
            { withCredentials: true }
          )
        ).subscribe({
          next: (res) => console.log(`[Backend] Loaded ${label}`, res),
          error: (err) => console.error(`[Backend] Failed to load ${label} ${path}`, err),
        });
      };

      if (paths.adataPath) loadPath('read_adata', paths.adataPath, 'adata');
      if (paths.genieFiltPath) loadPath('read_network_genie', paths.genieFiltPath, 'network_genie');
      if (paths.spongeFiltPath) loadPath('read_network_sponge', paths.spongeFiltPath, 'network_sponge');
      // for multiome data
      // if (paths.fragmentsFilePath && paths.fragmentsIndexPath) {
      //   loadPath('read_fragments', paths.fragmentsFilePath, 'fragments');
      //   loadPath('read_fragments_index', paths.fragmentsIndexPath, 'fragments_index');
      // }
    });
  }

  title = 'frontend';
  public formsOpen = false;

  public closeSidenavForms(): void {
    this.formsOpen = false;
  }

  public openForm(): void {
    this.formsOpen = true;
  }

  @HostListener('window:beforeunload', ['$event'])
  beforeUnloadHandler(_event: Event) {
    // Could add a notification that the current session is being deleted
    this.sessionService
      .callWithSession(() =>
        this.http.post(
          `${this.sessionService.apiUrl}/delete_session`,
          {},
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => console.log('[Backend] Deleted session', res),
        error: (err) =>
          console.error('[Backend] Failed to delete session', err),
      });
  }
}
