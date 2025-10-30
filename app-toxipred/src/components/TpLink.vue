<template>
  <a :href="href" class="tp-link" :class="[`tp-link--${underline}`, { 'tp-link--decolorized': inheritColor, 'tp-link--decolorized-hover': decolorizeHover }]">
    {{ label }}
  </a>
</template>

<script setup lang="ts">
interface Props {
  label: string;
  href: string;
  underline?: 'hover' | 'always' | 'none';
  inheritColor?: boolean;
  decolorizeHover?: boolean;
}

withDefaults(defineProps<Props>(), {
  underline: 'always'
});
</script>

<style scoped lang="scss">
.tp-link {
  color: var(--text-brand-regular);
  font-size: inherit;
  font-weight: inherit;
  text-decoration: underline;
  transition: all 0.2s ease-in-out;

  &:hover {
    color: var(--text-brand-medium-dark);
  }

  &:active {
    color: var(--text-brand-dark);
  }

  &--decolorized {
    color: inherit;
  }

  &--decolorized-hover {
    &:hover {
      color: inherit;
    }
  }

  &--always {
    text-decoration: underline;
    &:hover { text-decoration: underline; }
  }

  &--hover {
    text-decoration: none;
    &:hover { text-decoration: underline; }
  }

  &--none {
    text-decoration: none;
    &:hover { text-decoration: none; }
  }
}
</style>
