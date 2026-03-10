import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { HexagonPlotComponent } from './hexagon-plot/hexagon-plot.component';
import { InfoComponent } from './info/info.component';


export const routes: Routes = [
  { path: '', component: HexagonPlotComponent },
  { path: 'info', component: InfoComponent },
  { path: '**', redirectTo: '' }
];


//@NgModule({
//  imports: [RouterModule.forRoot(routes)],
//  exports: [RouterModule]
//})
//export class AppRoutingModule {}
