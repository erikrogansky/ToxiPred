<template>
  <tp-page class="row items-center justify-evenly">
    <div class="row items-center justify-between full-width">
      <h1 class="tp-main-heading">My Workspace</h1>
      <tp-button
        label="New prediction"
        @click="showNewPredictionDialog = true"
      />
    </div>

    <tp-advanced-filter
      :filters="filterConfig"
      :sort-options="sortConfig"
      :search-fields="searchFieldsConfig"
      default-sort-field="date"
      default-sort-direction="desc"
      default-search-field="all"
      @update:state="onFilterStateChange"
    />

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
        v-for="[id, j] in filteredAndSortedItems"
        :key="id"
        :id="id"
        :name="getDisplayName(j)"
        :prediction-type="getPredictionTypeLabel(j.model)"
        :result="getResultLabel(j.model, j.prediction)"
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

    <!-- New Prediction Dialog -->
    <tp-new-prediction-dialog v-model="showNewPredictionDialog" />
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpButton from 'components/TpButton.vue';
import TpPredictionRow from 'components/TpPredictionRow.vue';
import TpIcon from 'src/components/TpIcon.vue';
import TpNewPredictionDialog from 'src/components/TpNewPredictionDialog.vue';
import { TpAdvancedFilter } from 'src/components/filter';
import type { FilterConfig, SortOption, FilterOption, FilterState } from 'src/components/filter';

import { computed, ref, onMounted, onUnmounted } from 'vue'
import { useJobsStore } from 'src/stores/jobs-store';
import { useModelsStore, type TestType } from 'src/stores/models-store';
import { useRouter } from 'vue-router';

const router = useRouter()
const store = useJobsStore()
const modelsStore = useModelsStore()

// Dialog state
const showNewPredictionDialog = ref(false);

// Labels for display
const testTypeLabels: Record<TestType, string> = {
  'in_vitro': 'In Vitro',
  'in_vivo': 'In Vivo',
  'in_chemico': 'In Chemico',
}

// Filter configuration - dynamically built from model data
const filterConfig = computed<FilterConfig[]>(() => {
  // Build result options from model details
  const resultOptions = new Set<string>();
  Object.values(modelsStore.details).forEach(detail => {
    if (detail.positive_label) resultOptions.add(detail.positive_label);
    if (detail.negative_label) resultOptions.add(detail.negative_label);
  });

  return [
    {
      key: 'result',
      label: 'Result',
      options: Array.from(resultOptions).map(r => ({ label: r, value: r })),
    },
    {
      key: 'testType',
      label: 'Test Type',
      options: modelsStore.getTestTypes.map(t => ({
        label: testTypeLabels[t] || t,
        value: t,
      })),
    },
  ];
});

const sortConfig: SortOption[] = [
  { label: 'Date', value: 'date' },
  { label: 'Name', value: 'name' },
  { label: 'Result', value: 'result' },
];

// Search fields - simple categories to search within
const searchFieldsConfig: FilterOption[] = [
  { label: 'All', value: 'all' },
  { label: 'Name', value: 'name' },
  { label: 'SMILES', value: 'smiles' },
  { label: 'Formula', value: 'formula' },
  { label: 'Result', value: 'result' },
  { label: 'Method', value: 'method' },
];

// Filter state
const filterState = ref<FilterState>({
  filters: {},
  sort: 'date_desc',
  sortDirection: 'desc',
  searchQuery: '',
  searchField: 'all',
});

const onFilterStateChange = (state: FilterState) => {
  filterState.value = state;
};

const getPredictionTypeLabel = (modelName: string): string => {
  const detail = modelsStore.getModelDetail(modelName)
  if (detail?.test_type) {
    return testTypeLabels[detail.test_type] || detail.test_type
  }
  return modelName
}

const getResultLabel = (modelName: string, prediction: number | null): string => {
  const detail = modelsStore.getModelDetail(modelName)
  if (prediction === 1) {
    return detail?.positive_label || 'Positive'
  }
  return detail?.negative_label || 'Negative'
}

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

// Filtered and sorted items
const filteredAndSortedItems = computed(() => {
  let result = [...items.value];

  const { filters, sort, searchQuery, searchField } = filterState.value;

  // Apply filters
  const resultFilters = filters.result;
  if (resultFilters && resultFilters.length > 0) {
    result = result.filter(([, job]) => {
      const resultLabel = getResultLabel(job.model, job.prediction);
      return resultFilters.includes(resultLabel);
    });
  }

  const testTypeFilters = filters.testType;
  if (testTypeFilters && testTypeFilters.length > 0) {
    result = result.filter(([, job]) => {
      const detail = modelsStore.getModelDetail(job.model);
      return detail?.test_type && testTypeFilters.includes(detail.test_type);
    });
  }

  // Apply search
  if (searchQuery && searchQuery.trim() !== '') {
    const query = searchQuery.toLowerCase().trim();
    result = result.filter(([id, job]) => {
      const detail = modelsStore.getModelDetail(job.model);
      const resultLabel = getResultLabel(job.model, job.prediction);
      const testTypeLabel = detail?.test_type ? testTypeLabels[detail.test_type] : '';

      switch (searchField) {
        case 'all':
          return (
            job.name?.toLowerCase().includes(query) ||
            job.trivial_name?.toLowerCase().includes(query) ||
            job.canonical_smiles?.toLowerCase().includes(query) ||
            job.formula?.toLowerCase().includes(query) ||
            resultLabel.toLowerCase().includes(query) ||
            testTypeLabel.toLowerCase().includes(query) ||
            id.toLowerCase().includes(query)
          );
        case 'name':
          return (
            job.name?.toLowerCase().includes(query) ||
            job.trivial_name?.toLowerCase().includes(query)
          );
        case 'smiles':
          return job.canonical_smiles?.toLowerCase().includes(query);
        case 'formula':
          return job.formula?.toLowerCase().includes(query);
        case 'result':
          return resultLabel.toLowerCase().includes(query);
        case 'method':
          return testTypeLabel.toLowerCase().includes(query);
        default:
          return true;
      }
    });
  }

  // Apply sorting
  result.sort(([, a], [, b]) => {
    switch (sort) {
      case 'date_desc':
        return b.createdAt - a.createdAt;
      case 'date_asc':
        return a.createdAt - b.createdAt;
      case 'name_asc': {
        const nameA = (a.trivial_name || a.name || '').toLowerCase();
        const nameB = (b.trivial_name || b.name || '').toLowerCase();
        return nameA.localeCompare(nameB);
      }
      case 'name_desc': {
        const nameA = (a.trivial_name || a.name || '').toLowerCase();
        const nameB = (b.trivial_name || b.name || '').toLowerCase();
        return nameB.localeCompare(nameA);
      }
      case 'result_desc':
        return (b.prediction ?? -1) - (a.prediction ?? -1);
      case 'result_asc':
        return (a.prediction ?? -1) - (b.prediction ?? -1);
      default:
        return 0;
    }
  });

  return result;
});

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

const getDisplayName = (j: { trivial_name: string | null, canonical_smiles: string | null, name: string | null, formula: string | null }) => {
  let displayName: string;
  
  // If trivial_name is available, use it
  if (j.trivial_name) {
    displayName = j.trivial_name;
  } 
  // If no trivial_name but has canonical_smiles, use SMILES
  else if (j.canonical_smiles) {
    displayName = j.canonical_smiles;
  }
  // Fallback to name or 'Unknown compound'
  else {
    displayName = j.name || 'Unknown compound';
  }
  
  // Add formula if available
  if (j.formula) {
    displayName = displayName + ' | ' + formatChemFormula(j.formula);
  }
  
  return displayName;
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
