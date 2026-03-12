// This file contains environment-specific configuration values
// It is typically used for development mode in Angular projects.
// Angular replaces this file with `environment.prod.ts` during production builds
// (see `fileReplacements` array in `angular.json`).

export const environment = {
  production: false, // Indicates that this is NOT a production build (affects logging, debugging, etc.)

  // Base URL for API requests during development.
  // This is used by HttpClient calls in the app.
  apiBaseUrl: '',
};
