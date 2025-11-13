import { boot } from 'quasar/wrappers'

declare global {
  interface Window {
    RDKit?: unknown
    __rdkitReady?: Promise<unknown>
  }
}

export default boot(async () => {
  if (window.__rdkitReady) {
    await window.__rdkitReady
  }
  if (!window.RDKit) {
    console.error('RDKit not available after __rdkitReady.')
  } else {
    console.log('✅ RDKit ready in boot')
  }
})
