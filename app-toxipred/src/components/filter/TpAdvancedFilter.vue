<template>
  <div class="tp-advanced-filter">
    <!-- Left side: Filters and Sort -->
    <div class="tp-advanced-filter__left">
      <!-- Filters row -->
      <div v-if="filters.length > 0" class="tp-advanced-filter__row">
        <tp-icon icon-name="filter" weight="regular" :size="20" class="tp-advanced-filter__icon" />
        <div class="tp-advanced-filter__chips">
          <tp-filter-chip
            v-for="filter in filters"
            :key="filter.key"
            :label="filter.label"
            :options="filter.options"
            :model-value="filterValues[filter.key] || []"
            @update:model-value="(val) => updateFilter(filter.key, val)"
            @clear="clearFilter(filter.key)"
          />
        </div>
      </div>

      <!-- Sort row -->
      <div v-if="sortOptions.length > 0" class="tp-advanced-filter__row">
        <tp-icon icon-name="sort" weight="regular" :size="20" class="tp-advanced-filter__icon" />
        <tp-sort-control
          :options="sortOptions"
          :field="sortField"
          :direction="sortDirection"
          @update:field="sortField = $event"
          @update:direction="sortDirection = $event"
        />
      </div>
    </div>

    <!-- Right side: Search -->
    <div v-if="searchFields.length > 0" class="tp-advanced-filter__right">
      <tp-search-input
        :query="searchQuery"
        :field="searchField"
        :fields="searchFields"
        @update:query="searchQuery = $event"
        @update:field="searchField = $event"
      />
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, watch, reactive } from 'vue';
import TpIcon from 'src/components/TpIcon.vue';
import TpFilterChip from './TpFilterChip.vue';
import TpSortControl from './TpSortControl.vue';
import TpSearchInput from './TpSearchInput.vue';
import type { FilterConfig, SortOption, FilterOption, SortDirection, FilterState } from './types';
import { getSortValue } from './types';

const props = withDefaults(defineProps<{
  filters?: FilterConfig[];
  sortOptions?: SortOption[];
  searchFields?: FilterOption[];
  defaultSortField?: string;
  defaultSortDirection?: SortDirection;
  defaultSearchField?: string;
}>(), {
  filters: () => [],
  sortOptions: () => [],
  searchFields: () => [],
  defaultSortField: 'date',
  defaultSortDirection: 'desc',
  defaultSearchField: 'all',
});

const emit = defineEmits<{
  (e: 'update:state', state: FilterState): void;
}>();

// Reactive state
const filterValues = reactive<Record<string, string[]>>({});
const sortField = ref(props.defaultSortField);
const sortDirection = ref<SortDirection>(props.defaultSortDirection);
const searchQuery = ref('');
const searchField = ref(props.defaultSearchField);

// Filter actions
const updateFilter = (key: string, values: string[]) => {
  filterValues[key] = values;
};

const clearFilter = (key: string) => {
  filterValues[key] = [];
};

// Emit state on any change
const emitState = () => {
  emit('update:state', {
    filters: { ...filterValues },
    sort: getSortValue(sortField.value, sortDirection.value),
    sortDirection: sortDirection.value,
    searchQuery: searchQuery.value,
    searchField: searchField.value,
  });
};

// Watch all state changes
watch(
  [() => ({ ...filterValues }), sortField, sortDirection, searchQuery, searchField],
  emitState,
  { deep: true, immediate: true }
);
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
    gap: 8px;
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
    gap: 8px;
  }

  &__right {
    flex-shrink: 0;
  }
}
</style>
