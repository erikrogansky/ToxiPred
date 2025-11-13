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

    <div class="column full-width">
      <div
        class="row items-center justify-between full-width cursor-pointer q-py-sm"
        v-for="[id, j] in items"
        :key="id"
        @click="openJob(id)"
      >
        <!-- LEFT: name + formula -->
        <div class="column">
          <span>
            {{ j.name || 'Unnamed compound' }}
            <span v-if="j.formula">
              |
              <span v-html="formatChemFormula(j.formula)"></span>
            </span>
          </span>
          <span class="text-caption text-grey-7">
            {{ j.model }}
          </span>
        </div>

        <!-- MIDDLE: state + prediction -->
        <div class="column items-end q-mx-md">
          <span class="text-caption">
            {{ prettyState(j.state) }}
            <span v-if="j.percent !== null && !terminal(j.state)">
              · {{ j.percent }}%
            </span>
          </span>
          <span v-if="j.prediction !== null" class="text-caption">
            Prediction:
            <strong>{{ j.prediction === 1 ? 'Toxic' : 'Non-toxic' }}</strong>
          </span>
        </div>

        <!-- RIGHT: time + trash -->
        <div class="row items-center">
          <span class="text-caption q-mr-md">
            {{ formatTime(j.createdAt) }}
          </span>

          <tp-icon-button
            iconName="trash"
            weight="regular"
            :size="24"
            @click="remove(id)"
          />
        </div>
      </div>
    </div>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpButton from 'components/TpButton.vue';
import TpIcon from 'components/TpIcon.vue';
import TpIconButton from 'components/TpIconButton.vue';

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
const remove = (id: string) => store.remove(id)

const pendingPredictions = computed(() =>
  Array.from(store.jobs.values()).filter(j => !terminal(j.state))
)

const openJob = async (id: string) => {
  await router.push({ name: 'job-overview', params: { job_id: id } });
}

const goToNewPrediction = async () => {
  // adjust route name to whatever your home / input page is
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

const prettyState = (state: string) => {
  switch (state) {
    case 'PENDING': return 'Queued'
    case 'STARTED':
    case 'PROGRESS': return 'Running'
    case 'SUCCESS': return 'Finished'
    case 'FAILURE': return 'Failed'
    case 'REVOKED': return 'Cancelled'
    default: return state
  }
}
</script>

<style scoped lang="scss">
.tp-pending-panel {
  background-color: #FFF3E0;
  color: #EF6C00;
  padding: 8px 16px;
  border-radius: 8px;
  font-weight: 500;
  gap: 8px;
  margin-bottom: 24px;
}
</style>
