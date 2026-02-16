import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';

export interface Paths {
  adataMainPath?: string;
  adataComparePath?: string;
  genieFiltPath?: string;
  spongeFiltPath?: string;
  hexagonPath?: string;
}

@Injectable({
  providedIn: 'root'
})
export class PathsService {
  private paths = new BehaviorSubject<Paths>({});
  public paths$ = this.paths.asObservable();

  constructor() {}

  updatePaths(newPaths: Partial<Paths>): void {
    this.paths.next({
      ...this.paths.value,
      ...newPaths
    });
  }

  getPaths(): Paths {
    return this.paths.value;
  }
}
