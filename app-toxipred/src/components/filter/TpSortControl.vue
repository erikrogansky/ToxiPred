<template>
  <div class="tp-sort-control">
    <!-- Sort field chip -->
    <q-chip
      clickable
      outline
      dense
      :label="'Sort: ' + currentLabel"
      class="tp-sort-control__chip"
    >
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
    </q-chip>

    <!-- Direction buttons -->
    <q-chip
      clickable
      outline
      dense
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
@use 'src/css/helpers/glass' as *;

.tp-sort-control {
  display: flex;
  align-items: center;

  &__chip {
    @include glass-chip;
  }

  &__menu {
    @include glass-menu;
  }

  &__content {
    display: flex;
    flex-direction: column;
    min-width: 160px;
  }

  &__item {
    border-radius: 10px;
    padding: 8px 12px;
    transition: all 0.15s ease;

    &:hover {
      background-color: var(--surface-brand-extra-light);
      color: var(--text-brand-regular);
    }
  }

  &__direction {
    @include glass-chip;
    padding: 5px 8px !important;
    min-width: auto;

    &--active {
      @include glass-chip-active;
    }
  }
}
</style>
