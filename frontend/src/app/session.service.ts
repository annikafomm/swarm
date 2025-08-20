import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable, ReplaySubject, switchMap } from 'rxjs';

@Injectable({ providedIn: 'root' })
export class SessionService {
  // NOTE: If you access the website under 127.0.0.1 instead of
  // localhost you need to change this here as well.
  public apiUrl = 'http://localhost:3000';
  private sessionReady$ = new ReplaySubject<void>(1);

  constructor(private http: HttpClient) {}

  initSession(): void {
    const sessionId = crypto.randomUUID();
    this.http
      .post(
        `${this.apiUrl}/create_session/${sessionId}`,
        {},
        { withCredentials: true },
      )
      .subscribe({
        next: () => {
          console.log('[Backend] Created session');
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
