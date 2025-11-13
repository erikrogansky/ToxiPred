<template>
  <tp-page class="row items-center justify-evenly">
    <div class="row items-center justify-between full-width">
      <h1 class="tp-main-heading">My Workspace</h1>
      <tp-button
        label="New prediction"
        @click="goToNewPrediction"
      />
    </div>

    <div
      class="tp-pending-panel row items-center full-width"
      v-if="pendingPredictions.length > 0"
    >
      <tp-icon iconName="clock" weight="regular" :size="22" />
      <span>
        {{ pendingPredictions.length }}
        compound{{ pendingPredictions.length === 1 ? '' : 's' }}
        waiting in queue
      </span>
    </div>

    <q-list class="tp-prediction-list column items-center justify-between full-width">
      <tp-prediction-row 
        v-for="[id, j] in items"
        :key="id"
        :id="id"
        :name="getNameWithFormula(j.name, j.formula)"
        :result="j.prediction === 1 ? 'Toxic' : 'Non-toxic'"
        :time="formatTime(j.createdAt)"
        @row-click="openJob"                
        @selected-change="onSelectedChange"
        :selected="selectedIds.has(id)"
      />
    </q-list>

    <q-popup class="tp-selected-popup" v-if="selectedIds.size > 0">
      <span>Selected {{ selectedIds.size }} item{{ selectedIds.size > 1 ? 's' : '' }}</span>
      <tp-button size="small" variant="outline" label="Select all" @click="selectAll" />
      <tp-button size="small" variant="outline" label="Dselect all" @click="deselectAll" />
      <tp-button size="small" variant="outline" label="Delete selected" @click="deleteSelected" />
    </q-popup>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpButton from 'components/TpButton.vue';
import TpPredictionRow from 'components/TpPredictionRow.vue';

import { computed, ref, onMounted, onUnmounted } from 'vue'
import { useJobsStore } from 'src/stores/jobs-store';
import { useRouter } from 'vue-router';

const router = useRouter()
const store = useJobsStore()

const now = ref(Date.now())
let timer: number | null = null

onMounted(() => {
  now.value = Date.now()

  timer = window.setInterval(() => {
    now.value = Date.now()
  }, 60_000)
})

onUnmounted(() => {
  if (timer !== null) {
    window.clearInterval(timer)
  }
})

const items = computed(() => Array.from(store.jobs.entries()))
const terminal = (s: string) => ['SUCCESS', 'FAILURE', 'REVOKED'].includes(s)
const deleteSelected = () => {
  for (const id of selectedIds.value) {
    void store.remove(id)
    selectedIds.value.delete(id)
  }
}

const pendingPredictions = computed(() =>
  Array.from(store.jobs.values()).filter(j => !terminal(j.state))
)

const openJob = async (id: string) => {
  await router.push({ name: 'job-overview', params: { job_id: id } });
}

const goToNewPrediction = async () => {
  await router.push({ name: 'home' })
}

const subscriptMap: Record<string, string> = {
  '0': '₀',
  '1': '₁',
  '2': '₂',
  '3': '₃',
  '4': '₄',
  '5': '₅',
  '6': '₆',
  '7': '₇',
  '8': '₈',
  '9': '₉'
}

const getNameWithFormula = (name: string | null, formula: string | null) => {
  let final: string;
  final = name || 'Unknown compound';
  if (formula) {
    final = final + ' | ' + formatChemFormula(formula)
  }
  return final;
}

const formatChemFormula = (formula: string) => {
  if (!formula) return ''
  return formula.replace(/\d/g, digit => subscriptMap[digit] ?? digit)
}

const formatTime = (time: string | number) => {
  const parsed = new Date(time)
  if (Number.isNaN(parsed.getTime())) {
    return time
  }

  const diff = now.value - parsed.getTime()
  const absDiff = Math.abs(diff)

  const units = [
    { limit: 60_000, divisor: 1_000, unit: 'second' },
    { limit: 3_600_000, divisor: 60_000, unit: 'minute' },
    { limit: 86_400_000, divisor: 3_600_000, unit: 'hour' },
    { limit: 604_800_000, divisor: 86_400_000, unit: 'day' },
    { limit: 2_592_000_000, divisor: 604_800_000, unit: 'week' },
    { limit: 31_536_000_000, divisor: 2_592_000_000, unit: 'month' },
  ]

  for (const { limit, divisor, unit } of units) {
    if (absDiff < limit) {
      const value = Math.max(1, Math.floor(absDiff / divisor))
      const plural = value === 1 ? unit : `${unit}s`
      return diff >= 0 ? `${value} ${plural} ago` : `in ${value} ${plural}`
    }
  }

  const years = Math.max(1, Math.floor(absDiff / 31_536_000_000))
  const yearUnit = years === 1 ? 'year' : 'years'
  return diff >= 0 ? `${years} ${yearUnit} ago` : `in ${years} ${yearUnit}`
}


const selectedIds = ref(new Set<string>())

function onSelectedChange({ id, selected }: { id: string; selected: boolean }) {
  const next = new Set(selectedIds.value)
  if (selected) {
    next.add(id)
  } else {
    next.delete(id)
  }
  selectedIds.value = next
}

function selectAll() {
  selectedIds.value = new Set(items.value.map(([id]) => id))
}

function deselectAll() {
  selectedIds.value = new Set()
}

</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-pending-panel {
  background-color: #FFF3E0;
  color: #EF6C00;
  padding: 8px 16px;
  border-radius: 8px;
  font-weight: 500;
  gap: 8px;
  margin-bottom: 24px;
}

.tp-prediction-list {
  gap: 16px
}

.tp-selected-popup {
  position: fixed;
  padding: 16px;
  background: color-with-opacity(var(--surface-white), $opacity-high);
}
</style>
