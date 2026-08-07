import { Injectable } from '@angular/core';
import {
  HttpInterceptor,
  HttpRequest,
  HttpHandler,
  HttpEvent,
} from '@angular/common/http';
import { Observable } from 'rxjs';
import { environment } from '../environments/environment';

/**
 * Automatically prepends `environment.apiBaseUrl` to every relative URL
 * (i.e. URLs starting with '/') made through Angular's HttpClient.
 *
 * Dev:  apiBaseUrl = 'http://localhost:3000'
 *       '/api/datasets'  →  'http://localhost:3000/api/datasets'
 *
 * Prod: apiBaseUrl = '/swarm'
 *       '/api/datasets'  →  '/swarm/api/datasets'
 *       '/swarm/api/...' →  unchanged  (already prefixed)
 *
 * Absolute URLs (http:// / https://) are never modified.
 */
@Injectable()
export class ApiBaseUrlInterceptor implements HttpInterceptor {
  private readonly baseUrl = environment.apiBaseUrl;

  intercept(
    req: HttpRequest<unknown>,
    next: HttpHandler,
  ): Observable<HttpEvent<unknown>> {
    // Only transform relative URLs
    if (!this.baseUrl || !req.url.startsWith('/')) {
      return next.handle(req);
    }

    // Avoid double-prefixing when the URL already starts with the base path
    // (e.g. session.service.ts already prepends apiBaseUrl for some calls)
    const basePath = this.baseUrl.startsWith('/')
      ? this.baseUrl        // prod:  '/swarm'
      : undefined;          // dev:   'http://localhost:3000' — no relative prefix to skip

    if (basePath && req.url.startsWith(basePath)) {
      return next.handle(req);
    }

    const prefixedReq = req.clone({ url: `${this.baseUrl}${req.url}` });
    return next.handle(prefixedReq);
  }
}
