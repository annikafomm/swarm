import { Component } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import { HexagonPlotComponent } from "./hexagon-plot/hexagon-plot.component";

@Component({
  selector: 'app-root',
  imports: [RouterOutlet, HexagonPlotComponent],
  templateUrl: './app.component.html',
  styleUrl: './app.component.scss'
})
export class AppComponent {
  title = 'frontend';
}
