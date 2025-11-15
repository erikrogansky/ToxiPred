/* eslint-disable @typescript-eslint/prefer-promise-reject-errors */
import { boot } from 'quasar/wrappers'

type PlotlyLike = {
  newPlot: (
    el: HTMLElement,
    data: unknown[],
    layout?: Record<string, unknown>,
    config?: Record<string, unknown>,
  ) => unknown
  purge: (el: HTMLElement) => unknown
}

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
    Plotly?: PlotlyLike
    __plotlyReady?: Promise<void>
  }
}

export default boot(async () => {
  if (!window.__plotlyReady) {
    window.__plotlyReady = loadScript('https://cdn.plot.ly/plotly-2.32.0.min.js')
  }

  try {
    await window.__plotlyReady
  } catch (e) {
    console.error('Failed to load Plotly', e)
  }

  if (!window.Plotly) {
    console.error('Plotly failed to attach to window.Plotly')
  }
})
