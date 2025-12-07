<template>
  <div class="tp-sort-control">
    <!-- Sort field chip -->
    <q-chip
      clickable
      outline
      dense
      square
      :label="'Sort: ' + currentLabel"
      class="tp-sort-control__chip"
    />
    <q-menu class="tp-sort-control__menu">
      <div class="tp-sort-control__content">
        <q-item
          v-for="option in options"
          :key="option.value"
          clickable
          v-close-popup
          @click="$emit('update:field', option.value)"
          :active="field === option.value"
          class="tp-sort-control__item"
        >
          <q-item-section>{{ option.label }}</q-item-section>
        </q-item>
      </div>
    </q-menu>

    <!-- Direction buttons -->
    <q-chip
      clickable
      outline
      dense
      square
      class="tp-sort-control__direction"
      :class="{ 'tp-sort-control__direction--active': direction === 'asc' }"
      @click="$emit('update:direction', 'asc')"
    >
      <tp-icon icon-name="arrow-up-01" weight="regular" :size="16" />
    </q-chip>
    
    <q-chip
      clickable
      outline
      dense
      square
      class="tp-sort-control__direction"
      :class="{ 'tp-sort-control__direction--active': direction === 'desc' }"
      @click="$emit('update:direction', 'desc')"
    >
      <tp-icon icon-name="arrow-down-01" weight="regular" :size="16" />
    </q-chip>
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue';
import TpIcon from 'src/components/TpIcon.vue';
import type { SortOption, SortDirection } from './types';

const props = defineProps<{
  options: SortOption[];
  field: string;
  direction: SortDirection;
}>();

defineEmits<{
  (e: 'update:field', value: string): void;
  (e: 'update:direction', value: SortDirection): void;
}>();

const currentLabel = computed(() => {
  return props.options.find(o => o.value === props.field)?.label || 'Sort';
});
</script>

<style scoped lang="scss">
.tp-sort-control {
  display: flex;
  gap: 8px;
  align-items: center;

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
  }

  &__menu {
    background-color: var(--surface-white);
    border: 1px solid var(--stroke-extra-light);
    border-radius: 8px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
  }

  &__content {
    display: flex;
    flex-direction: column;
    padding: 8px;
    min-width: 150px;
  }

  &__item {
    border-radius: 4px;

    &:hover {
      background-color: var(--surface-brand-extra-light);
    }
  }

  &__direction {
    background-color: transparent;
    border: 1px solid var(--stroke-regular);
    color: var(--text);
    cursor: pointer;
    padding: 0 8px;
    min-width: auto;
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
}
</style>
