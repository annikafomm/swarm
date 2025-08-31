// paths.service.ts
import { Injectable } from '@angular/core';
import { BehaviorSubject } from 'rxjs';
import { DEFAULT_PATHS } from './constants';

@Injectable({
  providedIn: 'root'
})
export class PathsService {
  // BehaviorSubject so components can subscribe & update immediately
  private pathsSource = new BehaviorSubject(DEFAULT_PATHS);
  paths$ = this.pathsSource.asObservable();

  // get current snapshot
  get currentPaths() {
    return this.pathsSource.value;
  }

  // update paths (e.g. after upload)
  updatePaths(newPaths: typeof DEFAULT_PATHS) {
    this.pathsSource.next(newPaths);
  }
}
