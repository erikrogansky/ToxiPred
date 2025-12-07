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

  &__icon {
    color: var(--text-medium);
    flex-shrink: 0;
  }

  &__field {
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

  &__chip {
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

  &__menu {
    background-color: var(--surface-white);
    border: 1px solid var(--stroke-extra-light);
    border-radius: 8px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
  }

  &__menu-content {
    display: flex;
    flex-direction: column;
    min-width: 150px;
  }

  &__item {
    border-radius: 4px;

    &:hover {
      background-color: var(--surface-brand-extra-light);
    }
  }
}
</style>
