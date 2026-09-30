import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable, Subject, of } from 'rxjs';
import { catchError, map, shareReplay, tap } from 'rxjs/operators';
import { environment } from '../environments/environment';

/**
 * Ensembl gene ids: ENSG00000141510 (human), ENSMUSG00000059552 (mouse), optionally with a
 * version suffix. Kept in sync with `_ENSEMBL_ID_RE` in backend/main.py.
 */
const ENSEMBL_ID_PATTERN = /^ENS[A-Z]{0,4}G\d{6,}(?:\.\d+)?$/i;

/**
 * Translates Ensembl gene ids into the gene symbols a biologist actually reads.
 *
 * SPONGE is what makes this necessary. Its ceRNA modules, its interaction network and every
 * score derived from them are keyed by Ensembl id, so without translation the SPONGE half of
 * the UI shows `ENSG00000141510` in exactly the places the GENIE3 half shows `TP53` — module
 * dropdowns, network nodes, global and per-spot score tables alike.
 *
 * Rather than patching each of those call sites, everything that renders a gene identifier
 * goes through `label()` (or the `geneSymbol` pipe that wraps it). Identifiers that aren't
 * Ensembl ids — GENIE3's symbols, motif names, ligand-receptor pairs — pass through untouched,
 * so the helper is safe to apply to any gene-ish column without knowing where it came from.
 *
 * The map is deliberately **shared across datasets** rather than kept per dataset: a label is
 * a display concern with no dataset in scope at most call sites (a table cell, a d3 node), and
 * an Ensembl id means the same gene everywhere, so merging costs nothing and keeps `label()`
 * callable from anywhere. Where two datasets disagree — different Ensembl releases can rename
 * a gene — the first dataset loaded wins, and the raw id is still available via `idFor()`.
 */
@Injectable({ providedIn: 'root' })
export class GeneSymbolService {
  /** ENSG (upper-cased, both versioned and unversioned) -> symbol. */
  private readonly symbols = new Map<string, string>();

  /** In-flight / completed requests per dataset, so repeated selections don't refetch. */
  private readonly requests = new Map<string, Observable<Record<string, string>>>();

  private readonly changed = new Subject<void>();
  /**
   * Emits whenever newly loaded symbols became available. Template bindings pick the new
   * labels up on their own; imperative renderers (the d3 network graphs) subscribe to this to
   * redraw, since they'd otherwise keep whatever labels were current when they last ran.
   */
  public readonly changed$ = this.changed.asObservable();

  /**
   * Reads `environment.apiBaseUrl` directly rather than going through SessionService, which
   * exposes the same value. SessionService depends on DatasetService, and DatasetService is
   * what primes this service — taking the dependency would put this service inside that
   * import cycle for no gain, and a display-only helper has no business being entangled with
   * session bootstrapping.
   */
  private readonly apiUrl = environment.apiBaseUrl;

  constructor(private http: HttpClient) { }

  /**
   * Fetch and merge the symbol map for one dataset. Safe to call repeatedly — each dataset is
   * fetched at most once. A dataset with no Ensembl ids in its `var` yields an empty map, and
   * every id then simply keeps displaying as itself.
   */
  public load(datasetId?: string | null): Observable<Record<string, string>> {
    const cacheKey = datasetId || '__session__';
    const cached = this.requests.get(cacheKey);
    if (cached) {
      return cached;
    }

    const query = datasetId ? `?dataset_id=${encodeURIComponent(datasetId)}` : '';
    const request = this.http
      .get<{ symbols: Record<string, string> }>(
        `${this.apiUrl}/api/gene_symbols${query}`,
        { withCredentials: true },
      )
      .pipe(
        map((response) => response?.symbols ?? {}),
        tap((symbols) => this.merge(symbols)),
        catchError((err) => {
          // Not fatal: without a map the UI shows raw Ensembl ids, which is the old behaviour.
          console.warn(`[GeneSymbols] no symbol map for dataset "${cacheKey}"`, err);
          return of({} as Record<string, string>);
        }),
        shareReplay(1),
      );

    this.requests.set(cacheKey, request);
    // Kick the request off now rather than waiting for a subscriber: callers generally want the
    // map to be there by the time something renders, not to consume it directly.
    request.subscribe();
    return request;
  }

  private merge(symbols: Record<string, string>): void {
    let added = 0;
    for (const [id, symbol] of Object.entries(symbols || {})) {
      const key = String(id).trim().toUpperCase();
      if (!key || !symbol || this.symbols.has(key)) continue;
      this.symbols.set(key, String(symbol));
      added++;
    }
    if (added > 0) {
      this.changed.next();
    }
  }

  /** Whether `value` is an Ensembl gene id (and therefore worth translating). */
  public isEnsemblId(value: unknown): boolean {
    return typeof value === 'string' && ENSEMBL_ID_PATTERN.test(value.trim());
  }

  /**
   * The symbol for an Ensembl id, or null when the id is unknown — genes filtered out of the
   * dataset before it was written simply have no symbol to offer.
   */
  public symbolFor(value: unknown): string | null {
    if (!this.isEnsemblId(value)) return null;
    const key = String(value).trim().toUpperCase();
    return this.symbols.get(key) ?? this.symbols.get(key.split('.')[0]) ?? null;
  }

  /**
   * Display label for any gene identifier: the symbol when one is known, otherwise the value
   * unchanged. Non-Ensembl input is returned as-is, which is what makes this safe to apply
   * globally.
   */
  public label(value: unknown): string {
    if (typeof value !== 'string') return String(value ?? '');
    return this.symbolFor(value) ?? value;
  }

  /** `TP53 (ENSG00000141510)` — for the few places with room to show both. */
  public labelWithId(value: unknown): string {
    const symbol = this.symbolFor(value);
    return symbol ? `${symbol} (${String(value).trim()})` : this.label(value);
  }

  /**
   * The raw identifier behind a label, for tooltips and `title` attributes — so translating a
   * label never hides the id the underlying data is actually keyed by. Null when nothing was
   * substituted and the label already *is* the id.
   */
  public idFor(value: unknown): string | null {
    return this.symbolFor(value) ? String(value).trim() : null;
  }

  /** Whether `value` matches `query` by either its symbol or its raw id — for table filters. */
  public matches(value: unknown, query: string): boolean {
    const needle = query.trim().toLowerCase();
    if (!needle) return true;
    const raw = typeof value === 'string' ? value : String(value ?? '');
    if (raw.toLowerCase().includes(needle)) return true;
    const symbol = this.symbolFor(raw);
    return !!symbol && symbol.toLowerCase().includes(needle);
  }
}
