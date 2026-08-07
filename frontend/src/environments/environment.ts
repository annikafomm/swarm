// This file contains environment-specific configuration values
// It is typically used for development mode in Angular projects.
// Angular replaces this file with `environment.prod.ts` during production builds
// (see `fileReplacements` array in `angular.json`).

export const environment = {
  production: false,

  // In development (devcontainer), the backend runs directly on localhost:3000.
  // session.service.ts reads this value instead of hardcoding it.
  apiBaseUrl: 'http://localhost:3000',
};
