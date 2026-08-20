import { ApplicationConfig, provideZoneChangeDetection } from '@angular/core';
import { provideRouter, RouteReuseStrategy } from '@angular/router';
import { routes } from './app.routes';
import {
  provideHttpClient,
  withInterceptorsFromDi,
  HTTP_INTERCEPTORS,
} from '@angular/common/http';
import { DOCUMENT } from '@angular/common';
import { ComponentStateReuseStrategy } from './component-state-reuse.strategy';
import { ApiBaseUrlInterceptor } from './api-base-url.interceptor';

export const appConfig: ApplicationConfig = {
  providers: [
    provideZoneChangeDetection({ eventCoalescing: true }),
    provideRouter(routes),
    provideHttpClient(withInterceptorsFromDi()),
    { provide: DOCUMENT, useValue: document },
    { provide: RouteReuseStrategy, useClass: ComponentStateReuseStrategy },
    { provide: HTTP_INTERCEPTORS, useClass: ApiBaseUrlInterceptor, multi: true },
  ],
};

