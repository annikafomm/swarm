export const environment = {
  production: true,

  // Empty string = same-origin. In production, nginx serves the Angular app
  // and proxies /api/ requests to the backend container. No hardcoded host needed.
  apiBaseUrl: '',
};
