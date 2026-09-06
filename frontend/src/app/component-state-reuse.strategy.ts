import { Injectable } from '@angular/core';
import {
  RouteReuseStrategy,
  ActivatedRouteSnapshot,
  DetachedRouteHandle,
} from '@angular/router';

@Injectable()
export class ComponentStateReuseStrategy implements RouteReuseStrategy {
  private routeHandles = new Map<string, DetachedRouteHandle>();

  shouldDetach(route: ActivatedRouteSnapshot): boolean {
    // Only preserve the main visualization.
    return route.routeConfig?.path === '';
  }

  store(
    route: ActivatedRouteSnapshot,
    handle: DetachedRouteHandle | null
  ): void {
    if (handle) {
      this.routeHandles.set(this.getRoutePath(route), handle);
    }
  }

  shouldAttach(route: ActivatedRouteSnapshot): boolean {
    return this.routeHandles.has(this.getRoutePath(route));
  }

  retrieve(route: ActivatedRouteSnapshot): DetachedRouteHandle | null {
    return this.routeHandles.get(this.getRoutePath(route)) ?? null;
  }

  shouldReuseRoute(
    future: ActivatedRouteSnapshot,
    curr: ActivatedRouteSnapshot
  ): boolean {
    return future.routeConfig === curr.routeConfig;
  }

  private getRoutePath(route: ActivatedRouteSnapshot): string {
    return route.routeConfig?.path ?? '';
  }
}