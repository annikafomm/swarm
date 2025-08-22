import { Component } from '@angular/core';
import { Location } from '@angular/common';

@Component({
  selector: 'app-info',
  imports: [],
  templateUrl: './info.component.html',
  styleUrl: './info.component.scss'
})
export class InfoComponent {
  constructor(private location: Location) {}

  goBack() {
    this.location.back();
  }

}
