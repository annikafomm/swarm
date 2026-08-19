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
  // The server never exposes the raw session id (it lives in an httponly
  // cookie); `username` is the opaque, server-issued identity tied to it.
  public username: string | null = null;
  public isNewSession = false;

  constructor(private http: HttpClient, private datasetService: DatasetService) { }

  initSession(): void {
    this.http
      .post<{ username: string; is_new: boolean }>(
        `${this.apiUrl}/create_session`,
        {},
        { withCredentials: true },
      )
      .subscribe({
        next: (res) => {
          this.username = res.username;
          this.isNewSession = res.is_new;
          console.log(`[Backend] Session initialized (user: ${this.username}, new: ${this.isNewSession})`);
          // Load datasets after session is ready
          this.datasetService.loadAvailableDatasets();
          this.sessionReady$.next();
        },
        error: (err) =>
          console.error('[Backend] Failed to initialize session', err),
      });
  }

  callWithSession<T>(callback: () => Observable<T>): Observable<T> {
    return this.sessionReady$.pipe(switchMap(() => callback()));
  }
}
