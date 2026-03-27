import { Injectable } from '@angular/core';
import { BehaviorSubject } from 'rxjs';

/**
 * This service is required to access the geojson data from the download menu.
 * It may be a cleaner solution to load the data in the app component and input it
 * to the children, but this is simple and works for now.
 */
@Injectable({ providedIn: 'root' })
export class GeoDataService {
  private dataSource = new BehaviorSubject<any | null>(null);
  private dataSourceCompare = new BehaviorSubject<any | null>(null);
  data$ = this.dataSource.asObservable();
  dataCompare$ = this.dataSourceCompare.asObservable();

  setData(data: any) {
    this.dataSource.next(data);
  }

  setCompareData(data: any) {
    this.dataSourceCompare.next(data);
  }

  getData() {
    return this.dataSource.value;
  }

  getCompareData() {
    return this.dataSourceCompare.value;
  }
}
