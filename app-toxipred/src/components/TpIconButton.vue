<template>
  <button class="tp-icon-button" @click="$emit('click')">
      <i :class="`i-${iconName} ${weight}`" class="tp-icon-button__icon" :style="{ fontSize: `${size}px` }"></i>
  </button>
</template>

<script setup lang="ts">
import type { IconName, IconWeight } from 'src/css/icons/icon-names';

defineEmits(['click'])

defineProps<{
  iconName: IconName;
  weight: IconWeight<IconName>;
  size: number;
}>();

</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-icon-button{
  background: transparent;
  border: none;
  cursor: pointer;
  padding: 0;
  border-radius: 360px;
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 7px;

  &__icon {
    transition: all 0.3s ease-out;
  }

  &:hover {
    background-color: color-with-opacity(var(--surface-regular), $opacity-low);
  }

  &:active {
    background-color: color-with-opacity(var(--surface-regular), $opacity-medium);
  }
}

</style>
