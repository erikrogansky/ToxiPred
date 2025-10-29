<template>
  <button v-if="!href" class="tp-btn" :class="`tp-btn--${variant} tp-btn--${size}`" @click="$emit('click')">
    {{ label }}
  </button>
  <a v-else :href="href" class="tp-btn" :class="`tp-btn--${variant} tp-btn--${size}`">
    {{ label }}
  </a>
</template>

<script setup lang="ts">
interface Props {
  label: string;
  variant?: 'primary' | 'outline' | 'link';
  size?: 'small' | 'regular' | 'big';
  href?: string;
}

defineEmits(['click'])

withDefaults(defineProps<Props>(), {
  variant: 'primary',
  size: 'regular'
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-btn {
  font-weight: 900;
  border: none;
  border-radius: 360px;
  cursor: pointer;
  transition: background-color 0.3s ease;
  text-decoration: none;
  width: fit-content;

  &:focus-visible {
    outline: 2px solid var(--surface-brand-dark);
    outline-offset: 2px;
  }

  &--small {
    font-size: 12px;
    padding: 6px 12px;
  }

  &--regular {
    font-size: 14px;
    padding: 10px 20px;
  } 

  &--big {
    font-size: 16px;
     padding: 12px 24px;
  }

  &--primary {
    background-color: var(--surface-brand-regular);
    color: var(--text-inverse);

    &:hover {
      background-color: var(--surface-brand-medium-dark);
    }

    &:active {
      background-color: var(--surface-brand-dark);
    }
  }

  &--outline {
    background-color: transparent;
    color: var(--text);
    border: 1px solid var(--stroke-brand-regular);

    &:hover {
      background-color: color-with-opacity(var(--surface-brand-medium), $opacity-low);
    }

    &:active {
      background-color: color-with-opacity(var(--surface-brand-medium), $opacity-medium);
    }
  }

  &--link {
    background-color: transparent;
    color: var(--text);

    &:hover {
      background-color: color-with-opacity(var(--surface-regular), $opacity-low);
    }

    &:active {
      background-color: color-with-opacity(var(--surface-regular), $opacity-medium);
    }
  }
}
</style>
