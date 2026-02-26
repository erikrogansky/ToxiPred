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
.tp-sort-control {
  display: flex;
  align-items: center;
  gap: 4px;

  &__chip {
    background: var(--glass-background-light) !important;
    border: 1px solid var(--glass-border) !important;
    border-radius: 20px !important;
    color: var(--text) !important;
    cursor: pointer;
    transition: all 0.2s ease;
    font-weight: 500;
    font-size: 13px !important;
    padding: 5px 14px !important;
    height: 32px !important;
    -webkit-backdrop-filter: blur(var(--glass-blur));
    backdrop-filter: blur(var(--glass-blur));

    &:hover {
      border-color: var(--stroke-brand-regular) !important;
      background: var(--glass-background) !important;
      color: var(--text-brand-regular) !important;
      box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
    }
  }

  &__menu {
    background: var(--glass-background) !important;
    border: 1px solid var(--glass-border) !important;
    border-radius: 16px !important;
    box-shadow: var(--glass-shadow) !important;
    -webkit-backdrop-filter: blur(var(--glass-blur-strong)) !important;
    backdrop-filter: blur(var(--glass-blur-strong)) !important;
  }

  &__content {
    display: flex;
    flex-direction: column;
    padding: 8px;
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
    background: var(--glass-background-light) !important;
    border: 1px solid var(--glass-border) !important;
    border-radius: 20px !important;
    color: var(--text) !important;
    cursor: pointer;
    padding: 5px 8px !important;
    min-width: auto;
    height: 32px !important;
    transition: all 0.2s ease;
    -webkit-backdrop-filter: blur(var(--glass-blur));
    backdrop-filter: blur(var(--glass-blur));

    &:hover {
      border-color: var(--stroke-brand-regular) !important;
      background: var(--glass-background) !important;
      color: var(--text-brand-regular) !important;
      box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
    }

    &--active {
      border-color: var(--stroke-brand-regular) !important;
      background: var(--surface-brand-extra-light) !important;
      color: var(--text-brand-regular) !important;
      font-weight: 600;
    }
  }
}
</style>
