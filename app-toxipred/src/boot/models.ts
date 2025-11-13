import { boot } from 'quasar/wrappers'
import { useModelsStore } from 'src/stores/models-store'

export default boot(async () => {
  const store = useModelsStore()
  await store.init()
})
