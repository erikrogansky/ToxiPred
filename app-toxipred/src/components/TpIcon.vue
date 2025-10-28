<template>
  <i :class="`i-${iconName} ${weight}`" class="tp-icon" :style="{ fontSize: `${size}px` }"></i>
</template>

<script setup lang="ts">
import type { IconName, IconWeight } from 'src/css/icons/icon-names';

defineProps<{
  iconName: IconName;
  weight: IconWeight<IconName>;
  size: number;
}>();

</script>

<style scoped lang="scss">
.tp-icon {
  transition: all 0.3s ease-out;
}
</style>
