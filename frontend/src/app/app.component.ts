import { Component } from '@angular/core';
import { RouterOutlet, RouterLink } from '@angular/router';
import { CommonModule } from '@angular/common';
import { FormPageComponent } from './form-page/form-page.component';
//import { HexagonPlotComponent } from './hexagon-plot/hexagon-plot.component';

@Component({
  selector: 'app-root',
  standalone: true,
  imports: [CommonModule, FormPageComponent, RouterOutlet, RouterLink],
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
