import { Injectable } from '@angular/core';
import {
  RouteReuseStrategy,
  ActivatedRouteSnapshot,
  DetachedRouteHandle,
} from '@angular/router';

/**
 * Custom RouteReuseStrategy that preserves component state
 * when navigating between main page and info page.
 *
 * This ensures that when you navigate to /info and come back,
 * the HexagonPlotComponent maintains its state instead of being
 * destroyed and recreated.
 */
@Injectable()
export class ComponentStateReuseStrategy implements RouteReuseStrategy {
  private routeHandles = new Map<string, DetachedRouteHandle>();

  shouldDetach(route: ActivatedRouteSnapshot): boolean {
    // Always detach - we'll handle storing the component
    return true;
  }

  store(route: ActivatedRouteSnapshot, handle: DetachedRouteHandle | null): void {
    if (handle) {
      const routePath = this.getRoutePath(route);
      this.routeHandles.set(routePath, handle);
    }
  }

  shouldAttach(route: ActivatedRouteSnapshot): boolean {
    const routePath = this.getRoutePath(route);
    return this.routeHandles.has(routePath);
  }

  retrieve(route: ActivatedRouteSnapshot): DetachedRouteHandle | null {
    const routePath = this.getRoutePath(route);
    return this.routeHandles.get(routePath) || null;
  }

  shouldReuseRoute(
    future: ActivatedRouteSnapshot,
    curr: ActivatedRouteSnapshot
  ): boolean {
    // Reuse if the routes have the same path
    return future.routeConfig === curr.routeConfig;
  }

  private getRoutePath(route: ActivatedRouteSnapshot): string {
    return route.component ? (route.component as any).name : '';
  }
}
