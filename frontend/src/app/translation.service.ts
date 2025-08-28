import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { map, shareReplay } from 'rxjs/operators';
import { Observable } from 'rxjs';

@Injectable({
  providedIn: 'root',
})
export class TranslationService {
  private translations$!: Observable<Record<string, string>>;
  private translationsCache: Record<string, string> = {};

  constructor(private http: HttpClient) {
    this.translations$ = this.http
      .get<Record<string, string>>('assets/translations.json')
      .pipe(shareReplay(1));

    this.translations$.subscribe((dict) => {
      this.translationsCache = dict;
    });
  }

  getTranslations(): Observable<Record<string, string>> {
    return this.translations$;
  }

  translate(key: string): Observable<string> {
    return this.translations$.pipe(map((dict) => dict[key] ?? key));
  }

  translateSync(key: string): string {
    return this.translationsCache[key] ?? key;
  }
}
