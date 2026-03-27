import type { RouteRecordRaw } from 'vue-router';

const routes: RouteRecordRaw[] = [
  {
    path: '/',
    component: () => import('layouts/MainLayout.vue'),
    children: [{ path: '', component: () => import('pages/HomePage.vue'), name: 'home' }],
  },

  {
    path: '/workspace',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/WorkspacePage.vue'), name: 'workspace', },
      { path: 'job/:job_id', component: () => import('pages/JobOverviewPage.vue'), name: 'job-overview' },
    ],
  },

  {
    path: '/draw',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/KetcherPage.vue'), name: 'draw' },
    ],
  },

  {
    path: '/shared/:token',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/SharedJobPage.vue'), name: 'shared-job' },
    ],
  },

  {
    path: '/demos',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/DemosPage.vue'), name: 'demos' },
    ],
  },

  {
    path: '/privacy-policy',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/PrivacyPolicyPage.vue'), name: 'privacy-policy' },
    ],
  },

  {
    path: '/terms-of-service',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/TermsOfServicePage.vue'), name: 'terms-of-service' },
    ],
  },

  {
    path: '/cookie-policy',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/CookiePolicyPage.vue'), name: 'cookie-policy' },
    ],
  },

  {
    path: '/security',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/SecurityPage.vue'), name: 'security' },
    ],
  },

  {
    path: '/documentation',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', component: () => import('pages/DocumentationPage.vue'), name: 'documentation' },
    ],
  },

  {
    path: '/:catchAll(.*)*',
    component: () => import('layouts/MainLayout.vue'),
    children: [{ path: '', component: () => import('pages/ErrorNotFound.vue') }],
  },
];

export default routes;
