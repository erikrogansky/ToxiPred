<template>
  <div v-if="model" class="tp-nav__overlay" @click="close" />
  <nav class="tp-nav__panel" :class="{ 'is-open': model }" role="dialog" aria-modal="true" aria-label="Mobile navigation">
    <button class="tp-nav__close" aria-label="Close menu" @click="close">×</button>

    <ul class="tp-nav__list">
      <li v-for="(item, i) in items" :key="i" class="tp-nav__item">
        <a class="tp-nav__link" :href="item.href" @click="onItemClick()">{{ item.label }}</a>
      </li>
    </ul>
  </nav>
</template>

<script setup lang="ts">
import { onMounted, onBeforeUnmount, watch } from 'vue';

type NavItem = { label: string; href: string };

const model = defineModel<boolean>({ default: false });

const props = withDefaults(defineProps<{
  items?: NavItem[]
  closeOnNavigate?: boolean
}>(), {
  items: () => [],
  closeOnNavigate: true
});

function open()  { model.value = true; }
function close() { model.value = false; }
function toggle(){ model.value = !model.value; }
defineExpose({ open, close, toggle });

function onEsc(e: KeyboardEvent) {
  if (e.key === 'Escape' && model.value) close();
}
onMounted(() => document.addEventListener('keydown', onEsc));
onBeforeUnmount(() => document.removeEventListener('keydown', onEsc));

watch(model, (v) => {
  // lock scroll when open
  document.body.style.overflow = v ? 'hidden' : '';
});

function onItemClick() {
  if (props.closeOnNavigate) close();
}
</script>

<style scoped lang="scss">
.tp-nav__overlay {
  position: fixed;
  inset: 0;
  background: rgba(0,0,0,.5);
  z-index: 9998;
}

.tp-nav__panel {
  position: fixed;
  top: 0; right: 0; bottom: 0;
  width: min(82vw, 340px);
  background: var(--surface, #fff);
  box-shadow: -8px 0 24px rgba(0,0,0,.18);
  transform: translateX(100%);
  transition: transform .22s ease-out;
  z-index: 9999;
  display: flex;
  flex-direction: column;
  padding: 16px;
}

.tp-nav__panel.is-open {
  transform: translateX(0);
}

.tp-nav__close {
  border: 0;
  background: transparent;
  font-size: 28px;
  line-height: 1;
  align-self: flex-end;
  cursor: pointer;
}

.tp-nav__list {
  list-style: none;
  margin: 8px 0 0;
  padding: 0;
  display: grid;
  gap: 8px;
}

.tp-nav__item {}

.tp-nav__link {
  display: block;
  padding: 12px 8px;
  text-decoration: none;
  font-weight: 700;
  color: var(--text, #111);
}
.tp-nav__link:active { opacity: .8; }
</style>
