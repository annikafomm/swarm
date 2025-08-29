import { Component, OnInit, HostListener } from '@angular/core';
import { RouterOutlet, RouterLink } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
import { DownloadMenuComponent } from './download-menu/download-menu.component';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';

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
  ) {}

  ngOnInit() {
    this.sessionService.initSession();

    this.sessionService
      .callWithSession(() =>
        this.http.post(
          `${this.sessionService.apiUrl}/read_adata`,
          { path: '../backend/data/adata.h5ad' },
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => console.log('[Backend] Loaded adata', res),
        error: (err) => console.error('[Backend] Failed to load adata', err),
      });

    this.sessionService
      .callWithSession(() =>
        this.http.post(
          `${this.sessionService.apiUrl}/read_network_genie`,
          { path: '../backend/data/genie_network_filt.csv' },
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => console.log('[Backend] Loaded network_genie', res),
        error: (err) =>
          console.error('[Backend] Failed to load network_genie', err),
      });

    this.sessionService
      .callWithSession(() =>
        this.http.post(
          `${this.sessionService.apiUrl}/read_network_sponge`,
          { path: '../backend/data/sponge_network_filt.csv' },
          { withCredentials: true },
        ),
      )
      .subscribe({
        next: (res) => console.log('[Backend] Loaded network_sponge', res),
        error: (err) =>
          console.error('[Backend] Failed to load network_sponge', err),
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
