import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';

export interface Paths {
  adataPath?: string;
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
  private comparePaths = new BehaviorSubject<Paths>({});
  public comparePaths$ = this.comparePaths.asObservable();

  constructor() {}

  updatePaths(newPaths: Partial<Paths>, isCompare: boolean = false): void {
    if (isCompare) {
      this.comparePaths.next({
        ...this.comparePaths.value,
        ...newPaths
      });
    } else {
      this.paths.next({
        ...this.paths.value,
        ...newPaths
      });
    }
  }


  getPaths(isCompare: boolean = false): Paths {
    return isCompare ? this.comparePaths.value : this.paths.value;
  }
}
