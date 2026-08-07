import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable, ReplaySubject, switchMap } from 'rxjs';
import { DatasetService } from './datasets.service';
import { environment } from '../environments/environment';

@Injectable({ providedIn: 'root' })
export class SessionService {
  // Reads from environment.ts (dev: 'http://localhost:3000') or
  // environment.prod.ts (production: '' = same-origin, routed by nginx).
  // No more hardcoded URLs — no need for --skip-worktree on this file.
  public apiUrl = environment.apiBaseUrl;
  private sessionReady$ = new ReplaySubject<void>(1);
  public sessionId: string | null = null;


  constructor(private http: HttpClient, private datasetService: DatasetService) { }

  initSession(): void {
    const sessionId = crypto.randomUUID();
    this.sessionId = sessionId;
    this.http
      .post(
        `${this.apiUrl}/create_session/${sessionId}`,
        {},
        { withCredentials: true },
      )
      .subscribe({
        next: () => {
          console.log('[Backend] Created session');
          // Load datasets after session is ready
          this.datasetService.loadAvailableDatasets();
          this.sessionReady$.next();
        },
        error: (err) =>
          console.error('[Backend] Failed to create session', err),
      });


  }

  callWithSession<T>(callback: () => Observable<T>): Observable<T> {
    return this.sessionReady$.pipe(switchMap(() => callback()));
  }
}
