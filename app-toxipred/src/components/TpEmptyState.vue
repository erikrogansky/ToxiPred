<template>
  <div class="tp-empty-state" :class="{ 'tp-empty-state--compact': compact }">
    <div class="tp-empty-state__icon">
      <tp-icon :icon-name="icon" weight="regular" :size="compact ? 40 : 48" />
    </div>
    <h3 class="tp-empty-state__title">{{ title }}</h3>
    <p class="tp-empty-state__description">{{ description }}</p>
    <slot />
  </div>
</template>

<script setup lang="ts">
import TpIcon from './TpIcon.vue';
import type { IconName } from 'src/css/icons/icon-names';

withDefaults(defineProps<{
  icon?: IconName;
  title: string;
  description: string;
  compact?: boolean;
}>(), {
  icon: 'document',
  compact: false,
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/glass' as *;

.tp-empty-state {
  @include glass(var(--glass-background-light), var(--glass-blur), 20px);
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 12px;
  padding: 64px 32px;
  width: 100%;

  &--compact {
    padding: 48px 32px;
  }

  &__icon {
    color: var(--text-light);
    margin-bottom: 4px;
  }

  &__title {
    font-size: 18px;
    font-weight: 700;
    color: var(--text);
    margin: 0;
  }

  &__description {
    font-size: 14px;
    color: var(--text-medium);
    margin: 0 0 8px;
    text-align: center;
    max-width: 340px;
  }
}
</style>
