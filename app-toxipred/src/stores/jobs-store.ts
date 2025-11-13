import { defineStore } from 'pinia'
import { reactive } from 'vue'
import type { JobRecord, JobState } from 'src/types'
import { api } from 'src/boot/axios'

const LS_KEY = 'toxipred.jobs.v1'
const TERMINAL = new Set<JobState>(['SUCCESS', 'FAILURE', 'REVOKED'])

type JobSummary = Pick<JobRecord,
  'id' | 'model' | 'name' | 'formula' | 'state' | 'createdAt' | 'prediction'
>

export const useJobsStore = defineStore('jobs', () => {
  const jobs = reactive(new Map<string, JobRecord>())
  const streams = new Map<string, EventSource | null>()
  const pollers = new Map<string, number | null>()

  function toSummary(j: JobRecord): JobSummary {
    return {
      id: j.id,
      model: j.model,
      name: j.name,
      formula: j.formula,
      state: j.state,
      createdAt: j.createdAt,
      prediction: j.prediction,
    }
  }

  function save() {
    const arr = Array.from(jobs.values()).map(toSummary)
    localStorage.setItem(LS_KEY, JSON.stringify(arr))
  }

  function load() {
    const raw = localStorage.getItem(LS_KEY)
    if (!raw) return
    try {
      const arr: JobSummary[] = JSON.parse(raw)
      arr.forEach(s => {
        const record: JobRecord = {
          id: s.id,
          model: s.model,
          name: s.name ?? null,
          formula: s.formula ?? null,
          state: s.state,
          createdAt: s.createdAt,
          updatedAt: s.createdAt,

          percent: TERMINAL.has(s.state) ? 100 : null,
          msg: null,
          error: null,
          prediction: s.prediction ?? null,
        }
        jobs.set(s.id, record)
      })
    } catch {
      // ignore corrupted localStorage
    }
  }

  function upsert(j: Partial<JobRecord> & { id: string }) {
    const existing = jobs.get(j.id)
    const now = Date.now()

    const merged: JobRecord = {
      id: j.id,
      model: j.model ?? existing?.model ?? '',
      name: j.name ?? existing?.name ?? null,
      formula: j.formula ?? existing?.formula ?? null,
      state: j.state ?? existing?.state ?? 'PENDING',
      createdAt: existing?.createdAt ?? now,
      updatedAt: now,

      percent: j.percent ?? existing?.percent ?? null,
      msg: j.msg ?? existing?.msg ?? null,
      error: j.error ?? existing?.error ?? null,

      prediction: j.prediction ?? existing?.prediction ?? null,
    }

    jobs.set(j.id, merged)
    save()
  }

  async function remove(id: string) {
    stopTracking(id)
    jobs.delete(id)
    await api.delete(`/jobs/delete/${id}`)
    save()
  }

  function startTracking(id: string) {
    if (!('EventSource' in window)) return startPolling(id)
    stopTracking(id)

    const es = new EventSource(
      `${import.meta.env.API_URL || process.env.API_URL}/jobs/stream/${id}`,
    )
    streams.set(id, es)

    es.onmessage = (e) => {
      const p = JSON.parse(e.data)

      const res = p.result || null
      upsert({
        id,
        state: p.state,
        percent: p.progress_pct ?? (p.state === 'SUCCESS' ? 100 : null),
        msg: p.msg ?? null,
        error: p.error ?? null,
        name: res?.name ?? null,
        formula: res?.formula ?? null,
        model: res?.model ?? undefined,
        prediction: res?.prediction != null
          ? (Array.isArray(res.prediction) ? res.prediction[0] : res.prediction)
          : undefined,
      })

      if (TERMINAL.has(p.state)) stopTracking(id)
    }

    es.onerror = () => {
      es.close()
      streams.set(id, null)
      startPolling(id)
    }
  }

  function startPolling(id: string) {
    stopTracking(id)
    const timer = window.setInterval(() => {
      void (async () => {
        const res = await api.get(`/jobs/status/${id}`)
        const result = res.data.result || null

        upsert({
          id,
          state: res.data.state,
          percent: res.data.progress_pct ?? (res.data.state === 'SUCCESS' ? 100 : null),
          msg: res.data.msg ?? null,
          error: res.data.error ?? null,
          name: result?.name ?? null,
          formula: result?.formula ?? null,
          model: result?.model ?? undefined,
          prediction: result?.prediction != null
            ? (Array.isArray(result.prediction) ? result.prediction[0] : result.prediction)
            : undefined,
        })

        if (TERMINAL.has(res.data.state)) stopTracking(id)
      })()
    }, 900)
    pollers.set(id, timer)
  }

  function stopTracking(id: string) {
    const es = streams.get(id)
    if (es) es.close()
    streams.delete(id)

    const t = pollers.get(id)
    if (t) clearInterval(t)
    pollers.delete(id)
  }

  async function restoreAndAttach() {
    load()

    // 1) Reattach tracking for non-terminal jobs
    for (const j of jobs.values()) {
      if (!TERMINAL.has(j.state)) {
        startTracking(j.id)
      }
    }

    // 2) Ask backend which terminal jobs are no longer valid (missing/expired)
    const terminalIds = Array.from(jobs.values())
      .filter(j => TERMINAL.has(j.state))
      .map(j => j.id)

    if (terminalIds.length > 0) {
      try {
        const res = await api.post('/jobs/validate', {
          job_ids: terminalIds,
        })
        const invalidIds: string[] = res.data.invalid_ids ?? []

        invalidIds.forEach(id => {
          jobs.delete(id)
        })
      } catch (e) {
        console.error('Failed to validate jobs against backend:', e)
        // worst case: some stale entries stay in localStorage, no crash
      }
    }

    save()
  }

  async function fetchJobDetails(id: string) {
    const res = await api.get(`/jobs/result/${id}`)
    return res.data
  }

  return {
    jobs,
    upsert,
    remove,
    startTracking,
    restoreAndAttach,
    fetchJobDetails,
  }
})
