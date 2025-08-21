import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { HexagonPlotComponent } from './hexagon-plot/hexagon-plot.component';


export const routes: Routes = [
  { path: 'hexagon-plot', component: HexagonPlotComponent },
  { path: '', redirectTo: '/hexagon-plot', pathMatch: 'full' }
];


@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule {}
