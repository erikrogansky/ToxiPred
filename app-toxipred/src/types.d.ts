export type JobState = 'PENDING' | 'STARTED' | 'PROGRESS' | 'SUCCESS' | 'FAILURE' | 'REVOKED'

export interface JobRecord {
  id: string
  model: string
  name: string | null
  trivial_name: string | null
  formula: string | null
  state: JobState
  createdAt: number
  updatedAt: number

  percent: number | null
  msg: string | null
  error: string | null

  prediction: 0 | 1 | null
  canonical_smiles: string | null
  other_names: string[] | null
}

type JobSummary = Pick<JobRecord,
  'id' | 'model' | 'name' | 'trivial_name' | 'formula' | 'state' | 'createdAt' | 'prediction' | 'canonical_smiles'
>

declare module 'plotly.js-dist'
declare module 'plotly.js-dist-min'

// ketcher-standalone has types at dist/index.d.ts but package.json "exports" doesn't resolve them
declare module 'ketcher-standalone' {
  export class StandaloneStructServiceProvider {
    mode: string
    createStructService(options: Record<string, unknown>): unknown
  }
  export class StandaloneStructService {}
}

declare global {
  interface Window {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    $3Dmol?: any;
  }
}