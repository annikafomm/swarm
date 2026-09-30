import { Pipe, PipeTransform } from '@angular/core';
import { GeneSymbolService } from './gene-symbol.service';

/**
 * `{{ gene | geneSymbol }}` — renders an Ensembl gene id as its gene symbol, and anything else
 * unchanged. See GeneSymbolService for why this exists and why it is safe to apply broadly.
 *
 * Impure, like TranslatePipe: the map arrives asynchronously after the dataset is selected, and
 * an id already on screen has to pick up its symbol when it does. The work per call is a single
 * Map lookup.
 */
@Pipe({
  name: 'geneSymbol',
  standalone: true,
  pure: false,
})
export class GeneSymbolPipe implements PipeTransform {
  constructor(private geneSymbols: GeneSymbolService) { }

  transform(value: unknown): string {
    return this.geneSymbols.label(value);
  }
}
