<template>
  <div class="tp-search-input">
    <tp-icon 
      icon-name="search-normal" 
      weight="regular" 
      :size="18" 
      class="tp-search-input__icon" 
    />
    <input
      :value="query"
      @input="$emit('update:query', ($event.target as HTMLInputElement).value)"
      type="text"
      :placeholder="`Search in ${currentFieldLabel}...`"
      class="tp-search-input__field"
    />
    <q-chip
      clickable
      dense
      square
      :label="currentFieldLabel"
      class="tp-search-input__chip"
    >
      <q-menu class="tp-search-input__menu">
        <div class="tp-search-input__menu-content">
          <q-item
            v-for="option in fields"
            :key="option.value"
            clickable
            v-close-popup
            @click="$emit('update:field', option.value)"
            :active="field === option.value"
            class="tp-search-input__item"
          >
            <q-item-section>{{ option.label }}</q-item-section>
          </q-item>
        </div>
      </q-menu>
    </q-chip>
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue';
import TpIcon from 'src/components/TpIcon.vue';
import type { FilterOption } from './types';

const props = defineProps<{
  query: string;
  field: string;
  fields: FilterOption[];
}>();

defineEmits<{
  (e: 'update:query', value: string): void;
  (e: 'update:field', value: string): void;
}>();

const currentFieldLabel = computed(() => {
  return props.fields.find(f => f.value === props.field)?.label || 'All';
});
</script>

<style scoped lang="scss">
.tp-search-input {
  display: flex;
  align-items: center;
  gap: 10px;
  background: var(--glass-background-light);
  border: 1px solid var(--glass-border);
  border-radius: 20px;
  padding: 6px 14px;
  min-width: 300px;
  transition: all 0.2s ease;
  -webkit-backdrop-filter: blur(var(--glass-blur));
  backdrop-filter: blur(var(--glass-blur));

  &:hover {
    border-color: var(--stroke-brand-regular);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
  }

  &:focus-within {
    border-color: var(--stroke-brand-regular);
    background: var(--glass-background);
    box-shadow: 0 2px 12px rgba(0, 0, 0, 0.08);
  }

  &__icon {
    color: var(--text-brand-regular);
    flex-shrink: 0;
  }

  &__field {
    flex: 1;
    border: none;
    outline: none;
    background: transparent;
    font-size: 13px;
    font-weight: 500;
    color: var(--text);
    min-width: 0;

    &::placeholder {
      color: var(--text-light);
      font-weight: 400;
    }
  }

  &__chip {
    background: var(--surface-brand-extra-light);
    border: 1px solid var(--stroke-brand-regular);
    border-radius: 16px !important;
    color: var(--text-brand-regular);
    cursor: pointer;
    flex-shrink: 0;
    font-size: 11px;
    font-weight: 600;
    padding: 3px 10px !important;
    height: 24px !important;
    transition: all 0.2s ease;

    &:hover {
      background: var(--surface-brand-light);
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.06);
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

  &__menu-content {
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
}
</style>
