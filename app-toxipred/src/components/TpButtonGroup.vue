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
  border: 1px solid var(--glass-border);
  background: var(--glass-background-light);
  -webkit-backdrop-filter: blur(var(--glass-blur));
  backdrop-filter: blur(var(--glass-blur));
  width: 100%;

  &__indicator {
    position: absolute;
    top: -1px;
    bottom: -1px;
    left: 0;
    background: var(--glass-background);
    border: 1px solid var(--glass-border);
    border-radius: 360px;
    transition: transform 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    pointer-events: none;
    z-index: 0;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
    -webkit-backdrop-filter: blur(var(--glass-blur));
    backdrop-filter: blur(var(--glass-blur));
  }

  &__button {
    position: relative;
    z-index: 1;
    flex: 1;
    font-weight: 700;
    border: none;
    border-radius: 360px;
    cursor: pointer;
    transition: color 0.3s ease;
    text-decoration: none;
    background: transparent;
    padding: 8px 16px;
    font-size: 13px;
    color: var(--text-medium);

    &--active {
      color: var(--text);
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
