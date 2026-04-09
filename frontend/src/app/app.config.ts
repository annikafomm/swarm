import { ApplicationConfig, provideZoneChangeDetection } from '@angular/core';
import { provideRouter, RouteReuseStrategy } from '@angular/router';
import { routes } from './app.routes';
import {
  provideHttpClient,
  withInterceptorsFromDi,
} from '@angular/common/http';
import { ComponentStateReuseStrategy } from './component-state-reuse.strategy';

export const appConfig: ApplicationConfig = {
  providers: [
    provideZoneChangeDetection({ eventCoalescing: true }),
    provideRouter(routes),
    provideHttpClient(withInterceptorsFromDi()),
    { provide: RouteReuseStrategy, useClass: ComponentStateReuseStrategy },
  ],
};
