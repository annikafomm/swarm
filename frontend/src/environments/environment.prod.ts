export const environment = {
  production: true,

  // Sub-path deployment: global server nginx serves SWARM at /swarm/.
  // The global nginx strips /swarm/ before forwarding to the Docker nginx container,
  // so the inner nginx/FastAPI still receives clean paths (/api/..., /create_session/...).
  apiBaseUrl: '/swarm',
};
