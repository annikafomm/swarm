import { Component } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
import { HexagonPlotComponent } from './hexagon-plot/hexagon-plot.component';

@Component({
  selector: 'app-root',
  imports: [HexagonPlotComponent, CommonModule, FormPageComponent],
  templateUrl: './app.component.html',
  styleUrl: './app.component.scss'
})
export class AppComponent {
  title = 'frontend';
  public formsOpen = false;

   public closeSidenavForms(): void {
    this.formsOpen = false;
  }

  public openForm(): void {
    this.formsOpen = true;
  }
}
