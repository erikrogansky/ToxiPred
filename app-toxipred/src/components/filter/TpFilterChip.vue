<template>
  <q-chip
    clickable
    outline
    dense
    square
    :label="displayLabel"
    :removable="hasSelection"
    @remove="$emit('clear')"
    class="tp-filter-chip"
    :class="{ 'tp-filter-chip--active': hasSelection }"
  />
  <q-menu 
    class="tp-filter-chip__menu"
    @before-show="initPending"
  >
    <div class="tp-filter-chip__content">
      <q-checkbox
        v-for="option in options"
        :key="option.value"
        :label="option.label"
        :model-value="pendingValues.includes(option.value)"
        :true-value="true"
        :false-value="false"
        @update:model-value="(checked) => toggleOption(option.value, checked)"
        class="tp-filter-chip__checkbox"
      />
      <tp-button
        label="Apply"
        variant="primary"
        size="small"
        class="tp-filter-chip__apply"
        v-close-popup
        @click="apply"
      />
    </div>
  </q-menu>
</template>

<script setup lang="ts">
import { ref, computed } from 'vue';
import TpButton from 'src/components/TpButton.vue';
import type { FilterOption } from './types';

const props = defineProps<{
  label: string;
  options: FilterOption[];
  modelValue: string[];
}>();

const emit = defineEmits<{
  (e: 'update:modelValue', value: string[]): void;
  (e: 'clear'): void;
}>();

const pendingValues = ref<string[]>([]);

const hasSelection = computed(() => props.modelValue.length > 0);

const displayLabel = computed(() => {
  if (props.modelValue.length === 0) return props.label;
  const selectedLabels = props.options
    .filter(opt => props.modelValue.includes(opt.value))
    .map(opt => opt.label)
    .join(', ');
  return `${props.label}: ${selectedLabels}`;
});

const initPending = () => {
  pendingValues.value = [...props.modelValue];
};

const toggleOption = (value: string, checked: boolean) => {
  if (checked) {
    if (!pendingValues.value.includes(value)) {
      pendingValues.value.push(value);
    }
  } else {
    pendingValues.value = pendingValues.value.filter(v => v !== value);
  }
};

const apply = () => {
  emit('update:modelValue', [...pendingValues.value]);
};
</script>

<style scoped lang="scss">
.tp-filter-chip {
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

  &__menu {
    background-color: var(--surface-white);
    border: 1px solid var(--stroke-extra-light);
    border-radius: 8px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
  }

  &__content {
    display: flex;
    flex-direction: column;
    min-width: 180px;
  }

  &__checkbox {
    padding: 4px 8px;
  }

  &__apply {
    margin-top: 8px;
  }
}
</style>
