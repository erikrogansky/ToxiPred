export type JobState = 'PENDING' | 'STARTED' | 'PROGRESS' | 'SUCCESS' | 'FAILURE' | 'REVOKED'

export interface JobRecord {
  id: string
  model: string
  name: string | null
  formula: string | null
  state: JobState
  createdAt: number
  updatedAt: number

  percent: number | null
  msg: string | null
  error: string | null

  prediction: 0 | 1 | null
}

type JobSummary = Pick<JobRecord,
  'id' | 'model' | 'name' | 'formula' | 'state' | 'createdAt' | 'prediction'
>

declare module 'plotly.js-dist';
declare module 'plotly.js-dist-min';