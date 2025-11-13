/* eslint-disable @typescript-eslint/no-explicit-any */
/* eslint-disable @typescript-eslint/prefer-promise-reject-errors */
import { boot } from 'quasar/wrappers'

function loadScript(src: string) {
  return new Promise<void>((resolve, reject) => {
    const s = document.createElement('script')
    s.src = src
    s.async = true
    s.onload = () => resolve()
    s.onerror = (e) => reject(e)
    document.head.appendChild(s)
  })
}

declare global {
  interface Window {
    RDKit?: unknown
    __rdkitReady?: Promise<unknown>
    xsmiles?: any
    __xsmilesReady?: Promise<void>
  }
}

export default boot(async () => {
  if (window.__rdkitReady) {
    await window.__rdkitReady
  }
  if (!window.RDKit) {
    console.warn('XSMILES boot: RDKit still missing — XSMILES will load but won’t render until RDKit exists.')
  }

  if (!window.xsmiles) {
    await loadScript('/xsmiles/xsmiles.js')
  }
  if (!window.xsmiles) {
    console.error('XSMILES failed to attach to window.xsmiles')
    return
  }

  window.__xsmilesReady = Promise.resolve()
})
