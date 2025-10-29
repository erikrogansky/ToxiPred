<template>
  <div
    class="tp-dark-mode-toggle"
    :class="{ 'is-dark': darkMode }"
    role="switch"
    :aria-checked="darkMode"
    tabindex="0"
    @click="darkMode = !darkMode"
    @keydown.enter.prevent="darkMode = !darkMode"
    @keydown.space.prevent="darkMode = !darkMode"
  >
    <span class="knob" aria-hidden="true"></span>

    <span class="cell">
      <tp-icon icon-name="sun" weight="regular" :size="darkMode ? 19 : 20" :class="darkMode ? 'not-active' : 'active'" />
    </span>
    <span class="cell">
      <tp-icon icon-name="moon" weight="regular" :size="darkMode ? 20 : 19" :class="darkMode ? 'active' : 'not-active'" />
    </span>
  </div>
</template>

<script setup lang="ts">
import TpIcon from './TpIcon.vue'
import { ref, watch, onBeforeMount } from 'vue'

const darkMode = ref(false)

onBeforeMount(() => {
  const storedValue = localStorage.getItem('darkMode')
  if (storedValue) {
    darkMode.value = storedValue === 'true'
  } else {
    darkMode.value =
      window.matchMedia &&
      window.matchMedia('(prefers-color-scheme: dark)').matches
  }
  document.documentElement.setAttribute('data-theme', darkMode.value ? 'dark' : 'light')
})

watch(darkMode, (newValue) => {
  document.documentElement.setAttribute('data-theme', newValue ? 'dark' : 'light')
  localStorage.setItem('darkMode', String(newValue))
})
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-dark-mode-toggle {
  --cell: 28px;
  --pad: 4px;
  --gap: 6px;
  --ring: -1px;

  position: relative;
  display: inline-grid;
  grid-auto-flow: column;
  align-items: center;
  gap: var(--gap);
  padding: var(--pad);
  border: 1px solid color-with-opacity(var(--stroke-brand-regular), $opacity-regular);
  border-radius: 360px;
  cursor: pointer;
  user-select: none;
  -webkit-tap-highlight-color: transparent;
  will-change: transform;
  width: fit-content;

  &:focus-visible {
    outline: 2px solid var(--stroke-brand-regular);
    outline-offset: 2px;
  }

  .cell {
    width: var(--cell);
    height: var(--cell);
    display: grid;
    place-items: center;
    z-index: 1;

    .active {
      color: var(--text-dark);
    }

    .not-active {
      color: var(--text-brand-medium-dark);
    }
  }

  .knob {
    position: absolute;
    left: var(--ring);
    top: 50%;
    width: calc(var(--cell) + var(--pad) * 2 - var(--ring) * 2);
    height: calc(var(--cell) + var(--pad) * 2 - var(--ring) * 2);
    transform: translate(0, -50%);
    border-radius: 360px;
    background: var(--surface-brand-regular);
    box-shadow: 0 0 6px 1px rgba(23, 185, 136, 0.25);
    transition: transform 220ms cubic-bezier(.2,.7,.2,1);
  }

  &.is-dark .knob {
    transform: translate(calc(var(--cell) + var(--gap)), -50%);
  }

  @media (prefers-reduced-motion: reduce) {
    .knob { transition: none; }
  }
}
</style>
