<template>
  <div class="tp-descriptors-summary">
    <!-- Status Badges -->
    <div class="tp-status-badges">
      <div class="tp-badge tp-badge--result" :class="resultBadgeClass">
        <tp-icon 
          :icon-name="prediction === 0 ? 'tick-circle' : 'danger'" 
          :size="20"
          weight="regular"
        />
        <span>{{ resultLabel }}</span>
      </div>
      <div class="tp-badge tp-badge--confidence">
        <tp-icon 
          icon-name="percentage-circle" 
          :size="20"
          weight="regular"
        />
        <span>Confidence: {{ confidencePercent }}%</span>
      </div>
    </div>

    <!-- Descriptors Table -->
    <div class="tp-descriptors-table">
      <h3 class="tp-table-title">Descriptors summary</h3>
      
      <div class="tp-table-header">
        <div 
          class="tp-col tp-col--descriptor tp-col--sortable" 
          @click="toggleSort('descriptor')"
        >
          Descriptor
          <tp-icon 
            v-if="sortColumn === 'descriptor'" 
            icon-name="sort" 
            :size="14" 
            weight="regular" 
            :class="{ 'tp-sort-icon--desc': sortDirection === 'desc' }"
          />
        </div>
        <div 
          class="tp-col tp-col--type tp-col--sortable" 
          @click="toggleSort('type')"
        >
          Type
          <tp-icon 
            v-if="sortColumn === 'type'" 
            icon-name="sort" 
            :size="14" 
            weight="regular" 
            :class="{ 'tp-sort-icon--desc': sortDirection === 'desc' }"
          />
        </div>
        <div 
          class="tp-col tp-col--value tp-col--sortable" 
          @click="toggleSort('value')"
        >
          Value
          <tp-icon 
            v-if="sortColumn === 'value'" 
            icon-name="sort" 
            :size="14" 
            weight="regular" 
            :class="{ 'tp-sort-icon--desc': sortDirection === 'desc' }"
          />
        </div>
        <div 
          class="tp-col tp-col--importance tp-col--sortable" 
          @click="toggleSort('importance')"
        >
          Importance
          <tp-icon 
            v-if="sortColumn === 'importance'" 
            icon-name="sort" 
            :size="14" 
            weight="regular" 
            :class="{ 'tp-sort-icon--desc': sortDirection === 'desc' }"
          />
        </div>
      </div>

      <div class="tp-table-body">
        <div 
          v-for="descriptor in sortedDescriptors" 
          :key="descriptor.name"
          class="tp-table-row"
        >
          <div class="tp-col tp-col--descriptor">{{ descriptor.displayName }}</div>
          <div class="tp-col tp-col--type">
            <span class="tp-type-badge" :class="`tp-type-badge--${descriptor.typeClass}`">
              {{ descriptor.type }}
            </span>
          </div>
          <div class="tp-col tp-col--value">{{ descriptor.formattedValue }}</div>
          <div class="tp-col tp-col--importance">
            <div class="tp-stars">
              <span 
                v-for="star in 5" 
                :key="star" 
                class="tp-star"
                :class="{ 'tp-star--filled': star <= descriptor.stars }"
              >★</span>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed, ref } from 'vue';
import TpIcon from './TpIcon.vue';

interface Props {
  features: string[];
  values: (number | null)[];
  scores: number[];
  prediction: number | null;
  confidence: number | null;
  positiveLabel?: string;
  negativeLabel?: string;
}

const props = withDefaults(defineProps<Props>(), {
  positiveLabel: 'Toxic',
  negativeLabel: 'Non-toxic',
});

// Sorting state
type SortColumn = 'descriptor' | 'type' | 'value' | 'importance';
type SortDirection = 'asc' | 'desc';

const sortColumn = ref<SortColumn>('importance');
const sortDirection = ref<SortDirection>('desc');

function toggleSort(column: SortColumn) {
  if (sortColumn.value === column) {
    sortDirection.value = sortDirection.value === 'asc' ? 'desc' : 'asc';
  } else {
    sortColumn.value = column;
    sortDirection.value = column === 'importance' ? 'desc' : 'asc';
  }
}

// Descriptor type mapping
const descriptorTypes: Record<string, { type: string; typeClass: string; displayName: string }> = {
  // Physicochemical
  'qed': { type: 'Physicochemical', typeClass: 'physicochemical', displayName: 'QED' },
  'BalabanJ': { type: 'Physicochemical', typeClass: 'physicochemical', displayName: 'Balaban J' },
  'HallKierAlpha': { type: 'Physicochemical', typeClass: 'physicochemical', displayName: 'Hall-Kier Alpha' },
  
  // Electronic
  'MaxEStateIndex': { type: 'Electronic', typeClass: 'electronic', displayName: 'Max EState Index' },
  'MinAbsEStateIndex': { type: 'Electronic', typeClass: 'electronic', displayName: 'Min Abs EState Index' },
  'PEOE_VSA1': { type: 'Electronic', typeClass: 'electronic', displayName: 'PEOE VSA1' },
  'PEOE_VSA3': { type: 'Electronic', typeClass: 'electronic', displayName: 'PEOE VSA3' },
  'PEOE_VSA6': { type: 'Electronic', typeClass: 'electronic', displayName: 'PEOE VSA6' },
  'PEOE_VSA7': { type: 'Electronic', typeClass: 'electronic', displayName: 'PEOE VSA7' },
  'PEOE_VSA9': { type: 'Electronic', typeClass: 'electronic', displayName: 'PEOE VSA9' },
  'PEOE_VSA11': { type: 'Electronic', typeClass: 'electronic', displayName: 'PEOE VSA11' },
  'SMR_VSA6': { type: 'Electronic', typeClass: 'electronic', displayName: 'SMR VSA6' },
  'SlogP_VSA2': { type: 'Electronic', typeClass: 'electronic', displayName: 'SlogP VSA2' },
  'EState_VSA2': { type: 'Electronic', typeClass: 'electronic', displayName: 'EState VSA2' },
  'EState_VSA3': { type: 'Electronic', typeClass: 'electronic', displayName: 'EState VSA3' },
  'EState_VSA4': { type: 'Electronic', typeClass: 'electronic', displayName: 'EState VSA4' },
  'VSA_EState2': { type: 'Electronic', typeClass: 'electronic', displayName: 'VSA EState2' },
  'VSA_EState4': { type: 'Electronic', typeClass: 'electronic', displayName: 'VSA EState4' },
  'VSA_EState6': { type: 'Electronic', typeClass: 'electronic', displayName: 'VSA EState6' },
  
  // Topological
  'AvgIpc': { type: 'Topological', typeClass: 'topological', displayName: 'Avg IPC' },
  'BertzCT': { type: 'Topological', typeClass: 'topological', displayName: 'Bertz CT' },
  'Chi4n': { type: 'Topological', typeClass: 'topological', displayName: 'Chi4n' },
  
  // Structural
  'NumHeteroatoms': { type: 'Structural', typeClass: 'structural', displayName: 'Num Heteroatoms' },
  'fr_NH1': { type: 'Structural', typeClass: 'structural', displayName: 'NH1 Fragments' },
  'fr_NH2': { type: 'Structural', typeClass: 'structural', displayName: 'NH2 Fragments' },
  'fr_amide': { type: 'Structural', typeClass: 'structural', displayName: 'Amide Fragments' },
};

const resultLabel = computed(() => {
  if (props.prediction === null) return 'Unknown';
  return props.prediction === 0 ? props.negativeLabel : props.positiveLabel;
});

const resultBadgeClass = computed(() => {
  if (props.prediction === null) return '';
  return props.prediction === 0 ? 'tp-badge--negative' : 'tp-badge--positive';
});

const confidencePercent = computed(() => {
  if (props.confidence === null) return '—';
  return (props.confidence * 100).toFixed(1);
});

interface DescriptorRow {
  name: string;
  displayName: string;
  type: string;
  typeClass: string;
  value: number | null;
  formattedValue: string;
  score: number;
  stars: number;
}

const sortedDescriptors = computed<DescriptorRow[]>(() => {
  const len = Math.min(props.features.length, props.values.length, props.scores.length);
  const descriptors: DescriptorRow[] = [];

  // Find max absolute score for normalization
  let maxAbsScore = 0;
  for (let i = 0; i < len; i++) {
    const score = props.scores[i];
    if (score != null && !Number.isNaN(score)) {
      maxAbsScore = Math.max(maxAbsScore, Math.abs(score));
    }
  }

  for (let i = 0; i < len; i++) {
    const name = props.features[i] ?? '';
    const value = props.values[i] ?? null;
    const score = props.scores[i] ?? 0;

    const info = descriptorTypes[name] || {
      type: 'Other',
      typeClass: 'other',
      displayName: name,
    };

    // Calculate stars (1-5) based on relative importance
    const normalizedScore = maxAbsScore > 0 ? Math.abs(score) / maxAbsScore : 0;
    const stars = Math.max(1, Math.min(5, Math.ceil(normalizedScore * 5)));

    descriptors.push({
      name,
      displayName: info.displayName,
      type: info.type,
      typeClass: info.typeClass,
      value,
      formattedValue: formatValue(value, name),
      score,
      stars,
    });
  }

  // Sort based on current sort column and direction
  const dir = sortDirection.value === 'asc' ? 1 : -1;
  
  return descriptors.sort((a, b) => {
    switch (sortColumn.value) {
      case 'descriptor':
        return dir * a.displayName.localeCompare(b.displayName);
      case 'type':
        return dir * a.type.localeCompare(b.type);
      case 'value': {
        const aVal = a.value ?? 0;
        const bVal = b.value ?? 0;
        return dir * (aVal - bVal);
      }
      case 'importance':
      default:
        return dir * (Math.abs(a.score) - Math.abs(b.score));
    }
  });
});

function formatValue(value: number | null, name: string): string {
  if (value === null || Number.isNaN(value)) return '—';
  
  // Integer values for count-based descriptors
  if (name.startsWith('Num') || name.startsWith('fr_')) {
    return Math.round(value).toString();
  }
  
  // Small decimals
  if (Math.abs(value) < 0.01 && value !== 0) {
    return value.toExponential(2);
  }
  
  // Regular formatting
  if (Number.isInteger(value)) {
    return value.toString();
  }
  
  return value.toFixed(2);
}
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-descriptors-summary {
  overflow: hidden;
  width: 100%;
}

.tp-status-badges {
  display: flex;
  gap: 16px;
  margin-bottom: 24px;
}

.tp-badge {
  display: flex;
  align-items: center;
  gap: 8px;
  padding: 12px 24px;
  border-radius: 6px;
  font-weight: 600;
  font-size: 15px;
  flex: 1;
  justify-content: center;

  &--negative {
    background-color: #3FDBA5;
    color: #1a1a1a;
  }

  &--positive {
    background-color: #ef4444;
    color: white;
  }

  &--confidence {
    background-color: #F5C842;
    color: #1a1a1a;
  }
}

.tp-descriptors-table {
  background: color-with-opacity(var(--surface-brand-extra-light), $opacity-high);
  border-radius: 12px;
  padding: 20px;
}

.tp-table-title {
  font-size: 16px;
  font-weight: 600;
  margin-bottom: 16px;
  color: var(--text);
}

.tp-table-header {
  display: flex;
  //padding: 12px 16px;
  font-weight: 600;
  font-size: 13px;
  color: var(--text-brand-regular);
  background-color: rgba(128, 128, 128, 0.15);
  border-radius: 8px;
  margin-bottom: 4px;
}

.tp-col {
  padding: 12px 16px;
}

.tp-col--sortable {
  cursor: pointer;
  user-select: none;
  transition: background-color 0.15s ease;

  &:hover {
    background-color: rgba(128, 128, 128, 0.1);
    border-radius: 4px;
  }
}

.tp-table-body {
  max-height: 300px;
  overflow-y: auto;
}

.tp-table-row {
  display: flex;
  align-items: center;
  border-bottom: 1px solid rgba(0, 0, 0, 0.05);

  &:last-child {
    border-bottom: none;
  }

  &:hover {
    background-color: rgba(0, 0, 0, 0.02);
  }
}

.tp-col {
  padding: 12px 16px;
  text-align: left;

  &--descriptor {
    flex: 2;
    font-weight: 500;
    font-size: 14px;
  }

  &--type {
    flex: 1.5;
  }

  &--value {
    flex: 1;
    font-size: 14px;
  }

  &--importance {
    flex: 1.5;
    display: flex;
    align-items: center;
    gap: 4px;
  }
}

.tp-sort-icon--desc {
  transform: scaleY(-1);
}

.tp-type-badge {
  display: inline-block;
  padding: 4px 12px;
  border-radius: 50px;
  font-size: 12px;
  font-weight: 500;

  &--physicochemical {
    background-color: #E8B4F0;
    color: #5C2D6B;
  }

  &--electronic {
    background-color: #91D5F5;
    color: #1A5F7A;
  }

  &--structural {
    background-color: #A8E6CF;
    color: #2D5C4A;
  }

  &--topological {
    background-color: #FFBF87;
    color: #8B4513;
  }

  &--other {
    background-color: #D0D0D0;
    color: #4A4A4A;
  }
}

.tp-stars {
  display: flex;
  gap: 2px;
}

.tp-star {
  font-size: 16px;
  color: #D0D0D0;

  &--filled {
    color: #3FDBA5;
  }
}
</style>
