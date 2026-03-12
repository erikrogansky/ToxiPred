<template>
  <div class="tp-glass-card" :class="[`tp-glass-card--${variant}`, { 'tp-glass-card--hoverable': hoverable }]">
    <slot />
  </div>
</template>

<script setup lang="ts">
withDefaults(defineProps<{
  variant?: 'light' | 'regular' | 'elevated';
  hoverable?: boolean;
}>(), {
  variant: 'light',
  hoverable: false,
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/glass' as *;

.tp-glass-card {
  &--light {
    @include glass(var(--glass-background-light), var(--glass-blur), 14px);
  }

  &--regular {
    @include glass(var(--glass-background), var(--glass-blur), 14px);
  }

  &--elevated {
    @include glass-elevated(20px);
  }

  &--hoverable {
    transition: all 0.2s ease;

    &:hover {
      background: var(--glass-background);
      box-shadow: 0 4px 16px rgba(0, 0, 0, 0.06);
      border-color: var(--stroke-extra-light);
    }
  }
}
</style>
