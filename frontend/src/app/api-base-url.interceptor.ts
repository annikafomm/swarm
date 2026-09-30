import { Injectable } from '@angular/core';
import {
  HttpInterceptor,
  HttpRequest,
  HttpHandler,
  HttpEvent,
  HttpErrorResponse,
  HttpEventType,
} from '@angular/common/http';
import { Observable, throwError, EMPTY } from 'rxjs';
import { catchError, switchMap } from 'rxjs/operators';
import { environment } from '../environments/environment';

/**
 * Automatically prepends `environment.apiBaseUrl` to every relative URL
 * (i.e. URLs starting with '/') made through Angular's HttpClient.
 * Also automatically re-creates sessions on 403 "invalid session" errors (e.g. after backend restarts).
 */
@Injectable()
export class ApiBaseUrlInterceptor implements HttpInterceptor {
  private readonly baseUrl = environment.apiBaseUrl;

  intercept(
    req: HttpRequest<unknown>,
    next: HttpHandler,
  ): Observable<HttpEvent<unknown>> {
    let targetUrl = req.url;
    if (this.baseUrl && req.url.startsWith('/')) {
      const basePath = this.baseUrl.startsWith('/') ? this.baseUrl : undefined;
      if (!basePath || !req.url.startsWith(basePath)) {
        targetUrl = `${this.baseUrl}${req.url}`;
      }
    }

    const prefixedReq = req.clone({ url: targetUrl });

    return next.handle(prefixedReq).pipe(
      catchError((error: HttpErrorResponse) => {
        if (
          error.status === 403 &&
          error.error?.detail === 'invalid session' &&
          !targetUrl.includes('/create_session')
        ) {
          console.warn('[ApiBaseUrlInterceptor] Session invalidated by backend restart. Auto-creating fresh session...');
          const createSessionUrl = `${this.baseUrl}/create_session`;
          const createSessionReq = new HttpRequest('POST', createSessionUrl, {}, { withCredentials: true });

          return next.handle(createSessionReq).pipe(
            switchMap((event) => {
              if (event.type === HttpEventType.Response) {
                // Retry original request with newly issued cookie
                return next.handle(prefixedReq);
              }
              return EMPTY;
            })
          );
        }
        return throwError(() => error);
      })
    );
  }
}
