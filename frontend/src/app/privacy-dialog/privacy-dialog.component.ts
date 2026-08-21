import { Component, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MatDialogModule, MatDialogRef } from '@angular/material/dialog';
import { MatButtonModule } from '@angular/material/button';
import { MatIconModule } from '@angular/material/icon';
import { MatSlideToggleModule } from '@angular/material/slide-toggle';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-privacy-dialog',
  standalone: true,
  imports: [
    CommonModule,
    MatDialogModule,
    MatButtonModule,
    MatIconModule,
    MatSlideToggleModule,
    FormsModule,
  ],
  templateUrl: './privacy-dialog.component.html',
  styleUrl: './privacy-dialog.component.scss',
})
export class PrivacyDialogComponent implements OnInit {
  currentLang: 'en' | 'de' = 'en';
  isOptedOut: boolean = false;
  optOutStatusKnown: boolean = false;

  constructor(public dialogRef: MatDialogRef<PrivacyDialogComponent>) {}

  ngOnInit(): void {
    this.checkMatomoOptOutStatus();
  }

  setLanguage(lang: 'en' | 'de'): void {
    this.currentLang = lang;
  }

  checkMatomoOptOutStatus(): void {
    if (typeof window !== 'undefined' && (window as any)._paq) {
      (window as any)._paq.push([
        function (this: any) {
          const optedOut = typeof this.isUserOptedOut === 'function' ? this.isUserOptedOut() : false;
          setTimeout(() => {
            (window as any).__swarm_matomo_opted_out = optedOut;
          }, 0);
        },
      ]);

      setTimeout(() => {
        if ((window as any).__swarm_matomo_opted_out !== undefined) {
          this.isOptedOut = (window as any).__swarm_matomo_opted_out;
        } else {
          this.isOptedOut = localStorage.getItem('swarm_matomo_optout') === 'true';
        }
        this.optOutStatusKnown = true;
      }, 150);
    } else {
      this.isOptedOut = localStorage.getItem('swarm_matomo_optout') === 'true';
      this.optOutStatusKnown = true;
    }
  }

  onToggleOptOut(): void {
    if (typeof window !== 'undefined' && (window as any)._paq) {
      if (this.isOptedOut) {
        // User wants to opt out
        (window as any)._paq.push(['optUserOut']);
        localStorage.setItem('swarm_matomo_optout', 'true');
        (window as any).__swarm_matomo_opted_out = true;
      } else {
        // User wants to opt back in
        (window as any)._paq.push(['forgetUserOptOut']);
        localStorage.removeItem('swarm_matomo_optout');
        (window as any).__swarm_matomo_opted_out = false;
      }
    } else {
      if (this.isOptedOut) {
        localStorage.setItem('swarm_matomo_optout', 'true');
      } else {
        localStorage.removeItem('swarm_matomo_optout');
      }
    }
  }

  close(): void {
    this.dialogRef.close();
  }
}
