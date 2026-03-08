import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable, ReplaySubject, switchMap } from 'rxjs';
import { DatasetService } from './datasets.service';

@Injectable({ providedIn: 'root' })
export class SessionService {
  // NOTE: If you access the website under 127.0.0.1 instead of
  // localhost you need to change this here as well.
  public apiUrl = 'http://localhost:3005';
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
