// src/stores/models.ts
import { defineStore } from 'pinia'
import { api } from 'src/boot/axios'

export type TestType = 'in_vitro' | 'in_vivo' | 'in_chemico'
export type PredictionTarget = 'photo_irritation' | 'photo_toxicity' | 'corrosion'

export type ModelOption = {
  label: string
  value: string
}

const TEST_TYPE_LABELS: Record<TestType, string> = {
  in_vitro: 'In Vitro',
  in_vivo: 'In Vivo',
  in_chemico: 'In Chemico',
}

const PREDICTION_TARGET_LABELS: Record<PredictionTarget, string> = {
  photo_irritation: 'Photo Irritation',
  photo_toxicity: 'Phototoxicity',
  corrosion: 'Corrosion',
}

// Derive a short discriminator (e.g. "3T3", "3D", "in vivo") from the raw
// model name so dropdown labels stay unique when two models share the same
// test_type + prediction_target combination.
function modelDiscriminator(modelName: string): string | null {
  const paren = modelName.match(/\(([^)]+)\)\s*$/)
  if (paren && paren[1]) return paren[1]
  const parts = modelName.trim().split(/\s+/)
  const last = parts[parts.length - 1]
  if (!last) return null
  if (/^(XGB|GB|Corrosion|Phototox|Photo)$/i.test(last)) return null
  return last
}

export type ModelDetail = {
  file: string
  path: string
  kind: string
  features_in_spec?: number | null
  note?: string | null
  test_type?: TestType | null
  prediction_target?: PredictionTarget | null
  positive_label?: string | null
  negative_label?: string | null
  classification_threshold?: number | null
  dataset?: string | null
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
    
    // Get unique test types from available models
    getTestTypes: (s): TestType[] => {
      const types = new Set<TestType>()
      Object.values(s.details).forEach(detail => {
        if (detail.test_type) types.add(detail.test_type)
      })
      return Array.from(types)
    },
    
    // Get unique prediction targets from available models
    getPredictionTargets: (s): PredictionTarget[] => {
      const targets = new Set<PredictionTarget>()
      Object.values(s.details).forEach(detail => {
        if (detail.prediction_target) targets.add(detail.prediction_target)
      })
      return Array.from(targets)
    },
    
    // Get prediction targets available for a specific test type
    getPredictionTargetsForTestType: (s) => (testType: TestType): PredictionTarget[] => {
      const targets = new Set<PredictionTarget>()
      Object.values(s.details).forEach(detail => {
        if (detail.test_type === testType && detail.prediction_target) {
          targets.add(detail.prediction_target)
        }
      })
      return Array.from(targets)
    },
    
    // Find model by test type and prediction target
    getModelBySelection: (s) => (testType: TestType, predictionTarget: PredictionTarget): string | null => {
      for (const [modelName, detail] of Object.entries(s.details)) {
        if (detail.test_type === testType && detail.prediction_target === predictionTarget) {
          return modelName
        }
      }
      return null
    },

    // Formatted dropdown options with dedupe (appends a discriminator when
    // two models share the same test_type + prediction_target combo).
    getModelOptions: (s): ModelOption[] => {
      const counts = new Map<string, number>()
      const base = s.models.map((modelName) => {
        const detail = s.details[modelName]
        const parts: string[] = []
        if (detail?.prediction_target) {
          parts.push(PREDICTION_TARGET_LABELS[detail.prediction_target] || detail.prediction_target)
        }
        if (detail?.test_type) {
          parts.push(TEST_TYPE_LABELS[detail.test_type] || detail.test_type)
        }
        const label = parts.length > 0 ? parts.join(' - ') : modelName
        counts.set(label, (counts.get(label) ?? 0) + 1)
        return { modelName, label }
      })

      return base.map(({ modelName, label }) => {
        let finalLabel = label
        if ((counts.get(label) ?? 0) > 1) {
          const disc = modelDiscriminator(modelName)
          if (disc) finalLabel = `${label} (${disc})`
        }
        return { label: finalLabel, value: modelName }
      })
    },
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
