<template>
  <q-chip
    clickable
    outline
    dense
    :label="displayLabel"
    class="tp-filter-chip"
    :class="{ 'tp-filter-chip--active': hasSelection }"
  >
    <div v-if="hasSelection" class="tp-filter-chip__remove" @click.stop="$emit('clear')">
      <tp-icon icon-name="close" weight="regular" :size="14" />
    </div>
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
  </q-chip>
</template>

<script setup lang="ts">
import { ref, computed } from 'vue';
import TpButton from 'src/components/TpButton.vue';
import TpIcon from 'src/components/TpIcon.vue';
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

  &--active {
    border-color: var(--stroke-brand-regular) !important;
    background: var(--surface-brand-extra-light) !important;
    color: var(--text-brand-regular) !important;
    font-weight: 600;
  }

  &__menu {
    background: var(--glass-background) !important;
    border: 1px solid var(--glass-border) !important;
    border-radius: 16px !important;
    box-shadow: var(--glass-shadow) !important;
    -webkit-backdrop-filter: blur(var(--glass-blur-strong)) !important;
    backdrop-filter: blur(var(--glass-blur-strong)) !important;
    padding: 8px !important;
  }

  &__content {
    display: flex;
    flex-direction: column;
    min-width: 200px;
    gap: 4px;
  }

  &__checkbox {
    padding: 8px 12px;
    border-radius: 10px;
    transition: background-color 0.15s ease;

    &:hover {
      background-color: var(--surface-brand-extra-light);
    }
  }

  &__apply {
    margin: 8px;
    width: calc(100% - 16px);
  }

  &__remove {
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    margin-left: 4px;
    transition: opacity 0.2s ease;

    &:hover {
      opacity: 0.7;
    }
  }
}
</style>
