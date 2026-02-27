<template>
  <div class="tp-advanced-filter">
    <!-- Left side: Filters and Sort -->
    <div class="tp-advanced-filter__left">
      <!-- First row: Filter icon + filter chips -->
      <div class="tp-advanced-filter__row">
        <tp-icon icon-name="filter" weight="regular" :size="20" class="tp-advanced-filter__icon" />
        <div class="tp-advanced-filter__chips">
          <div v-for="filter in filters" :key="filter.key">
            <q-chip
              clickable
              outline
              dense
              square
              :label="getFilterLabel(filter)"
              :removable="hasActiveFilter(filter.key)"
              @remove="clearFilter(filter.key)"
              class="tp-advanced-filter__chip"
              :class="{ 'tp-advanced-filter__chip--active': hasActiveFilter(filter.key) }"
            />
            <q-menu 
              class="tp-advanced-filter__menu"
              @before-show="initPendingFilter(filter.key)"
            >
              <div class="tp-advanced-filter__menu-content">
                <q-checkbox
                  v-for="option in filter.options"
                  :key="option.value"
                  :label="option.label"
                  :model-value="pendingFilters[filter.key]?.includes(option.value) ?? false"
                  :true-value="true"
                  :false-value="false"
                  @update:model-value="(val) => togglePendingFilterOption(filter.key, option.value, val)"
                  class="tp-advanced-filter__checkbox"
                />
                <tp-button
                  label="Apply filter"
                  variant="primary"
                  size="small"
                  class="tp-advanced-filter__apply-btn"
                  v-close-popup
                  @click="applyFilter(filter.key)"
                />
              </div>
            </q-menu>
          </div>
        </div>
      </div>

      <!-- Second row: Sort icon + sort chip + direction chip -->
      <div class="tp-advanced-filter__row">
        <tp-icon icon-name="sort" weight="regular" :size="20" class="tp-advanced-filter__icon" />
        <div class="tp-advanced-filter__chips">
          <q-chip
            clickable
            outline
            dense
            square
            :label="getSortLabel()"
            class="tp-advanced-filter__chip"
          />
          <q-menu class="tp-advanced-filter__menu">
            <div class="tp-advanced-filter__menu-content">
              <q-item
                v-for="option in sortOptions"
                :key="option.value"
                clickable
                v-close-popup
                @click="selectedSortField = option.value"
                :active="selectedSortField === option.value"
                class="tp-advanced-filter__sort-item"
              >
                <q-item-section>{{ option.label }}</q-item-section>
              </q-item>
            </div>
          </q-menu>
          
          <q-chip
            clickable
            outline
            dense
            square
            class="tp-advanced-filter__chip tp-advanced-filter__direction-chip"
            :class="{ 'tp-advanced-filter__chip--active': sortDirection === 'asc' }"
            @click="sortDirection = 'asc'"
          >
            <tp-icon 
              icon-name="arrow-up-01" 
              weight="regular" 
              :size="16" 
            />
          </q-chip>
          
          <q-chip
            clickable
            outline
            dense
            square
            class="tp-advanced-filter__chip tp-advanced-filter__direction-chip"
            :class="{ 'tp-advanced-filter__chip--active': sortDirection === 'desc' }"
            @click="sortDirection = 'desc'"
          >
            <tp-icon 
              icon-name="arrow-down-01" 
              weight="regular" 
              :size="16" 
            />
          </q-chip>
        </div>
      </div>
    </div>

    <!-- Right side: Search input with search-in chip -->
    <div class="tp-advanced-filter__right">
      <div class="tp-advanced-filter__search-wrapper">
        <tp-icon icon-name="search-normal" weight="regular" :size="18" class="tp-advanced-filter__search-icon" />
        <input
          v-model="searchQuery"
          type="text"
          :placeholder="`Search in ${getSearchInLabel()}...`"
          class="tp-advanced-filter__search-input"
        />
        <q-chip
          clickable
          dense
          square
          :label="getSearchInLabel()"
          class="tp-advanced-filter__search-chip"
        >
          <q-menu class="tp-advanced-filter__menu">
            <div class="tp-advanced-filter__menu-content">
              <q-item
                v-for="option in searchInOptions"
                :key="option.value"
                clickable
                v-close-popup
                @click="selectedSearchIn = option.value"
                :active="selectedSearchIn === option.value"
                class="tp-advanced-filter__sort-item"
              >
                <q-item-section>{{ option.label }}</q-item-section>
              </q-item>
            </div>
          </q-menu>
        </q-chip>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import TpButton from 'src/components/TpButton.vue';
import TpIcon from 'src/components/TpIcon.vue';

import { ref, reactive, watch, computed } from 'vue';

export interface FilterOption {
  label: string;
  value: string;
}

export interface FilterConfig {
  key: string;
  label: string;
  options: FilterOption[];
}

export interface SortOption {
  label: string;
  value: string;
}

export interface SearchInOption {
  label: string;
  value: string;
}

export interface FilterState {
  filters: Record<string, string[]>;
  sort: string;
  searchQuery: string;
  searchIn: string;
}

const props = withDefaults(defineProps<{
  filters?: FilterConfig[];
  sortOptions?: SortOption[];
  searchInOptions?: SearchInOption[];
  defaultSort?: string;
}>(), {
  filters: () => [],
  sortOptions: () => [
    { label: 'Date (Newest first)', value: 'date_desc' },
    { label: 'Date (Oldest first)', value: 'date_asc' },
  ],
  searchInOptions: () => [
    { label: 'All', value: 'all' },
  ],
  defaultSort: 'date_desc',
});

const emit = defineEmits<{
  (e: 'update:filterState', state: FilterState): void;
}>();

// State - applied filters
const selectValues = ref<Record<string, string[]>>({});
// State - pending filters (before applying)
const pendingFilters = reactive<Record<string, string[]>>({});

const selectedSortField = ref<string>(props.sortOptions[0]?.value || 'date');
const sortDirection = ref<'asc' | 'desc'>('desc');
const selectedSearchIn = ref<string>(props.searchInOptions[0]?.value || 'all');
const searchQuery = ref<string>('');

// Computed sort value combining field and direction
const selectedSort = computed(() => `${selectedSortField.value}_${sortDirection.value}`);



// Emit filter state whenever anything changes
const emitFilterState = () => {
  emit('update:filterState', {
    filters: { ...selectValues.value },
    sort: selectedSort.value,
    searchQuery: searchQuery.value,
    searchIn: selectedSearchIn.value,
  });
};

// Watch for changes and emit
watch([selectValues, selectedSortField, sortDirection, searchQuery, selectedSearchIn], emitFilterState, { deep: true });

// Initialize pending filter when menu opens
const initPendingFilter = (filterKey: string) => {
  pendingFilters[filterKey] = [...(selectValues.value[filterKey] || [])];
};

// Toggle pending filter option
const togglePendingFilterOption = (filterKey: string, optionValue: string, isChecked: boolean) => {
  if (!pendingFilters[filterKey]) {
    pendingFilters[filterKey] = [];
  }
  
  if (isChecked) {
    if (!pendingFilters[filterKey].includes(optionValue)) {
      pendingFilters[filterKey].push(optionValue);
    }
  } else {
    pendingFilters[filterKey] = pendingFilters[filterKey].filter(v => v !== optionValue);
  }
};

// Apply filter from pending to actual values
const applyFilter = (filterKey: string) => {
  selectValues.value[filterKey] = [...(pendingFilters[filterKey] || [])];
};

// Clear a specific filter
const clearFilter = (filterKey: string) => {
  selectValues.value[filterKey] = [];
  pendingFilters[filterKey] = [];
};

// Check if filter has active selections
const hasActiveFilter = (filterKey: string): boolean => {
  return (selectValues.value[filterKey]?.length ?? 0) > 0;
};

// Helper functions
const getFilterLabel = (filter: FilterConfig): string => {
  const selected = selectValues.value[filter.key];
  if (selected && selected.length > 0) {
    const selectedLabels = filter.options
      ?.filter(opt => selected.includes(opt.value))
      .map(opt => opt.label)
      .join(', ');
    return `${filter.label}: ${selectedLabels}`;
  }
  return filter.label;
};

const getSortLabel = (): string => {
  const option = props.sortOptions.find(opt => opt.value === selectedSortField.value);
  return 'Sort: ' + option?.label || 'Sort';
};

const getSearchInLabel = (): string => {
  const option = props.searchInOptions.find(opt => opt.value === selectedSearchIn.value);
  return `In: ${option?.label || 'All'}`;
};
</script>

<style scoped lang="scss">
.tp-advanced-filter {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  width: 100%;
  gap: 24px;

  &__left {
    display: flex;
    flex-direction: column;
  }

  &__row {
    display: flex;
    align-items: center;
    gap: 12px;
  }

  &__icon {
    color: var(--text);
    flex-shrink: 0;
  }

  &__chips {
    display: flex;
    flex-wrap: wrap;
  }

  &__chip {
    background-color: transparent;
    border: 1px solid var(--stroke-regular);
    color: var(--text);
    cursor: pointer;
    transition: all 0.2s ease;

    &:hover {
      border-color: var(--stroke-brand-regular);
      color: var(--text-brand-regular);
    }

    &--active {
      border-color: var(--stroke-brand-regular);
      background-color: var(--surface-brand-extra-light);
      color: var(--text-brand-regular);
    }
  }

  &__direction-chip {
    padding: 0 8px;
    min-width: auto;
  }

  &__menu {
    background-color: var(--surface-white);
    border: 1px solid var(--stroke-extra-light);
    border-radius: 8px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
  }

  &__menu-content {
    display: flex;
    flex-direction: column;
    min-width: 180px;
  }

  &__checkbox {
    padding: 4px 8px;
  }

  &__apply-btn {
    margin-top: 8px;
    align-self: stretch;
  }

  &__sort-item {
    border-radius: 4px;

    &:hover {
      background-color: var(--surface-brand-extra-light);
    }
  }

  &__right {
    flex-shrink: 0;
  }

  &__search-wrapper {
    display: flex;
    align-items: center;
    gap: 8px;
    background-color: var(--surface-white);
    border: 1px solid var(--stroke-regular);
    border-radius: 8px;
    padding: 6px 12px;
    min-width: 280px;
    transition: border-color 0.2s ease;

    &:focus-within {
      border-color: var(--stroke-brand-regular);
    }
  }

  &__search-icon {
    color: var(--text-medium);
    flex-shrink: 0;
  }

  &__search-input {
    flex: 1;
    border: none;
    outline: none;
    background: transparent;
    font-size: 14px;
    color: var(--text);
    min-width: 0;

    &::placeholder {
      color: var(--text-medium);
    }
  }

  &__search-chip {
    background-color: var(--surface-brand-extra-light);
    border: 1px solid var(--stroke-brand-regular);
    color: var(--text-brand-regular);
    cursor: pointer;
    flex-shrink: 0;
    font-size: 12px;
    transition: all 0.2s ease;

    &:hover {
      background-color: var(--surface-brand-light);
    }
  }
}
</style>