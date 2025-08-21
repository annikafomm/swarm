import { Component, OnInit, HostListener } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
import { HexagonPlotComponent } from './hexagon-plot/hexagon-plot.component';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';

@Component({
  selector: 'app-root',
  imports: [HexagonPlotComponent, CommonModule, FormPageComponent],
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
