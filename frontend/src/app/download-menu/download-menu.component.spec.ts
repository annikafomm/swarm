import { ComponentFixture, TestBed } from '@angular/core/testing';

import { DownloadMenuComponent } from './download-menu.component';

describe('DownloadMenuComponent', () => {
  let component: DownloadMenuComponent;
  let fixture: ComponentFixture<DownloadMenuComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [DownloadMenuComponent]
    })
    .compileComponents();

    fixture = TestBed.createComponent(DownloadMenuComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
