import { ComponentFixture, TestBed } from '@angular/core/testing';

import { HexagonPlotComponent } from './hexagon-plot.component';

describe('HexagonPlotComponent', () => {
  let component: HexagonPlotComponent;
  let fixture: ComponentFixture<HexagonPlotComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [HexagonPlotComponent]
    })
    .compileComponents();

    fixture = TestBed.createComponent(HexagonPlotComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
