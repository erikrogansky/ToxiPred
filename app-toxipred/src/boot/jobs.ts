// src/boot/jobs.ts
import { boot } from 'quasar/wrappers'
import { useJobsStore } from 'src/stores/jobs-store'

export default boot(async () => {
  const store = useJobsStore()
  await store.restoreAndAttach()
})
