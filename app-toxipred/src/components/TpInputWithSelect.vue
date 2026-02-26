<template>
  <div class="tp-input-with-select-wrapper">
    <div class="tp-input-with-select">
      <input
        v-model="inputValue"
        @keyup.enter="$emit('enter')"
        type="text"
        :placeholder="label"
        class="tp-input-with-select__field"
      />
      <div class="tp-input-with-select__chip">
        <span class="tp-input-with-select__chip-text">{{ currentLabel }}</span>
        <tp-icon icon-name="arrow-down-02" weight="regular" :size="14" class="tp-input-with-select__chip-icon" />
        <q-menu class="tp-input-with-select__menu">
          <div class="tp-input-with-select__menu-content">
            <q-item
              v-for="option in options"
              :key="option.value"
              clickable
              v-close-popup
              @click="selectedValue = option.value"
              :active="selectedValue === option.value"
              class="tp-input-with-select__item"
            >
              <q-item-section>
                <q-item-label class="paragraph-small bold">{{ option.label }}</q-item-label>
              </q-item-section>
            </q-item>
          </div>
        </q-menu>
      </div>
    </div>
    <div class="tp-input-with-select__hint" v-if="hint">{{ hint }}</div>
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue';
import TpIcon from 'src/components/TpIcon.vue';

const inputValue = defineModel<string>('inputValue')
const selectedValue = defineModel<string>('selectedValue')

const props = defineProps<{
  options: Array<{ label: string; value: string }>;
  hint: string;
  label: string;
  selectLabel: string;
}>();

defineEmits<{
  (e: 'enter'): void;
}>();

const currentLabel = computed(() => {
  const found = props.options.find(o => o.value === selectedValue.value);
  return found?.label || props.selectLabel;
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-input-with-select-wrapper {
  max-width: 700px;
  width: 100%;
}

.tp-input-with-select {
  display: flex;
  align-items: center;
  gap: 10px;
  border: 1px solid var(--glass-border);
  border-radius: 28px;
  padding: 10px 16px;
  width: 100%;
  background: var(--glass-background-light);
  -webkit-backdrop-filter: blur(var(--glass-blur));
  backdrop-filter: blur(var(--glass-blur));
  box-shadow: var(--glass-shadow);
  transition: all 0.25s ease;

  &:hover {
    border-color: var(--stroke-light);
    box-shadow: var(--glass-shadow-elevated);
  }

  &:focus-within {
    border-color: var(--stroke-brand-regular);
    box-shadow: 0 0 0 3px rgba(var(--color-brand-rgb, 0, 128, 96), 0.12), var(--glass-shadow-elevated);
  }

  &__field {
    flex: 1;
    border: none;
    outline: none;
    background: transparent;
    font-size: 14px;
    color: var(--text);
    min-width: 150px;

    &::placeholder {
      color: var(--text-medium);
    }
  }

  &__chip {
    background: var(--surface-brand-extra-light);
    border: 1px solid var(--stroke-brand-regular);
    color: var(--text-brand-regular);
    cursor: pointer;
    flex-shrink: 0;
    font-size: 13px;
    font-weight: 600;
    padding: 5px 12px;
    border-radius: 18px;
    transition: all 0.2s ease;
    display: flex;
    align-items: center;
    gap: 6px;
    user-select: none;

    &:hover {
      background: var(--surface-brand-light);
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.06);
    }
  }

  &__chip-text {
    max-width: 180px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }

  &__chip-icon {
    flex-shrink: 0;
    margin-left: 2px;
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
    min-width: 200px;
    padding: 6px;
  }

  &__item {
    border-radius: 10px;
    padding: 10px 12px;
    transition: all 0.15s ease;

    &:hover {
      background: var(--surface-brand-extra-light);
    }

    &.q-item--active {
      background: var(--surface-brand-extra-light);
      color: var(--text-brand-regular);
      font-weight: 600;
    }
  }
}

.tp-input-with-select__hint {
  font-size: 12px;
  color: var(--text-medium);
  padding-left: 16px;
  margin-top: 4px;
}
</style>
