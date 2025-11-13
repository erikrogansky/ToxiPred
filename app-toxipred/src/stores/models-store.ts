// src/stores/models.ts
import { defineStore } from 'pinia'
import { api } from 'src/boot/axios'

export type ModelDetail = {
  file: string
  path: string
  kind: string
  features_in_spec?: number | null
  note?: string | null
}

export type ModelsResponse = {
  available_models: string[]
  details: Record<string, ModelDetail>
}

type ModelsState = {
  models: string[]
  details: Record<string, ModelDetail>
  etag?: string | undefined
  lastFetched?: number | undefined
  isLoading: boolean
  _initialized: boolean
}

const CACHE_KEY = 'models-cache-v1'
const STALE_MS = 60 * 60 * 1000

export const useModelsStore = defineStore('models', {
  state: (): ModelsState => ({
    models: [],
    details: {},
    etag: undefined,
    lastFetched: undefined,
    isLoading: false,
    _initialized: false,
  }),

  getters: {
    getModels: (s) => s.models,
    getDetails: (s) => s.details,
    getModelDetail: (s) => (name: string) => s.details[name],
    isStale: (s) => !s.lastFetched || (Date.now() - s.lastFetched) > STALE_MS,
  },

  actions: {
    setModels(payload: ModelsResponse, etag?: string) {
      this.models = payload.available_models ?? []
      this.details = payload.details ?? {}
      this.lastFetched = Date.now()
      if (etag) this.etag = etag

      // localStorage.setItem(
      //   CACHE_KEY,
      //   JSON.stringify({
      //     models: this.models,
      //     details: this.details,
      //     etag: this.etag,
      //     lastFetched: this.lastFetched,
      //   })
      // )
    },

    loadCache() {
      try {
        const raw = localStorage.getItem(CACHE_KEY)
        if (!raw) return
        const { models, details, etag, lastFetched } = JSON.parse(raw)
        this.models = models ?? []
        this.details = details ?? {}
        this.etag = etag
        this.lastFetched = lastFetched
      } catch {
        /* ignore bad cache */
      }
    },

    async init() {
      if (this._initialized) return
      this._initialized = true

      this.loadCache()
      await this.revalidate()

      window.addEventListener('visibilitychange', () => {
        if (document.visibilityState === 'visible' && this.isStale) {
          void this.revalidate()
        }
      })
    },

    async revalidate() {
      if (this.isLoading) return
      this.isLoading = true
      try {
        const headers: Record<string, string> = {}
        if (this.etag) headers['If-None-Match'] = this.etag

        const resp = await api.get<ModelsResponse>('/models', {
          headers,
          validateStatus: () => true,
        })

        console.log(resp.data)

        if (resp.status === 304) {
          this.lastFetched = Date.now()
          return
        }

        if (resp.status >= 200 && resp.status < 300 && resp.data) {
          const etag = resp.headers['etag'] as string | undefined
          this.setModels(resp.data, etag)
        }
      } finally {
        this.isLoading = false
      }
    },

    async invalidateAndReload() {
      this.lastFetched = 0
      await this.revalidate()
    },
  },
})
