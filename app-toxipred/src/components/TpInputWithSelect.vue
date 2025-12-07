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
  gap: 8px;
  border: 1px solid rgba(0, 0, 0, 0.24);
  border-radius: 28px;
  padding: 8px 14px;
  width: 100%;
  transition: border-color 0.2s ease;

  &:focus-within {
    border-color: var(--stroke-brand-regular);
  }

  &:hover {
    border-color: rgba(0, 0, 0, 0.87);
  }

  &:focus-within:hover {
    border-color: var(--stroke-brand-regular);
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
    background-color: var(--surface-brand-extra-light);
    border: 1px solid var(--stroke-brand-regular);
    color: var(--text-brand-regular);
    cursor: pointer;
    flex-shrink: 0;
    font-size: 13px;
    padding: 4px 10px;
    border-radius: 16px;
    transition: all 0.2s ease;
    display: flex;
    align-items: center;
    gap: 4px;
    user-select: none;

    &:hover {
      background-color: var(--surface-brand-light);
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
    background-color: var(--surface-white);
    border: 1px solid var(--stroke-extra-light);
    border-radius: 12px;
    box-shadow: 0 4px 16px rgba(0, 0, 0, 0.12);
  }

  &__menu-content {
    display: flex;
    flex-direction: column;
    min-width: 200px;
  }

  &__item {
    border-radius: 8px;
    padding: 10px 12px;

    &:hover {
      background-color: var(--surface-brand-extra-light);
    }

    &.q-item--active {
      background-color: var(--surface-brand-light);
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
