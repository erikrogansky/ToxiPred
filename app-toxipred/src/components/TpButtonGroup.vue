<template>
  <div class="tp-btn-group">
    <div 
      class="tp-btn-group__indicator"
      :style="indicatorStyle"
    ></div>
    <button
      v-for="(label, index) in labels"
      :key="index"
      class="tp-btn-group__button"
      :class="{ 'tp-btn-group__button--active': index === activeIndex }"
      @click="$emit('click', index)"
      :disabled="disabled"
    >
      {{ label }}
    </button>
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue';

interface Props {
  labels: string[];
  activeIndex?: number;
  disabled?: boolean;
}

defineEmits(['click'])

const props = withDefaults(defineProps<Props>(), {
  activeIndex: 0,
  disabled: false
});

const indicatorStyle = computed(() => {
  const width = 100 / props.labels.length;
  const translateAmount = props.activeIndex * 100;
  
  return {
    transform: `translateX(${translateAmount}%)`,
    width: `${width}%`
  };
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-btn-group {
  position: relative;
  display: flex;
  border-radius: 1000px;
  border: 1px solid var(--stroke-extra-light);
  background: transparent;
  width: 100%;

  &__indicator {
    position: absolute;
    top: 0;
    bottom: 0;
    left: 0;
    background: white;
    border-radius: 360px;
    transition: transform 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    pointer-events: none;
    z-index: 0;
  }

  &__button {
    position: relative;
    z-index: 1;
    flex: 1;
    font-weight: 900;
    border: none;
    border-radius: 360px;
    cursor: pointer;
    transition: color 0.3s ease;
    text-decoration: none;
    background: transparent;
    padding: 8px 16px;
    color: var(--text-secondary, #666);

    &--active {
      color: var(--text-primary, #000);
    }

    &:disabled {
      cursor: not-allowed;
      opacity: 0.5;
    }

    &:focus-visible {
      outline: 2px solid var(--surface-brand-dark);
      outline-offset: 2px;
    }

    &:hover:not(:disabled) {
      color: var(--text-primary, #000);
    }
  }
}
</style>
