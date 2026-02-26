<template>
    <q-item
      v-ripple
      clickable
      @click="onRowClick"
      class="tp-prediction-row row items-center justify-between full-width g"
    >
      <q-item-section avatar>
        <q-checkbox
          :model-value="selected"
          @click.stop
          @update:model-value="val => emit('selected-change', { id: props.id, selected: val })"
        />
      </q-item-section>

      <q-item-section>
        <q-item-label>{{ name }}</q-item-label>
      </q-item-section>

      <q-item-section class="tp-prediction-row__prediction-type">
        <q-item-label>{{ predictionType }}</q-item-label>
      </q-item-section>

      <q-item-section class="tp-prediction-row__result">
        <q-item-label>{{ result }}</q-item-label>
      </q-item-section>

      <q-item-section class="tp-prediction-row__time">
        <q-item-label>{{ time }}</q-item-label>
      </q-item-section>

      <q-item-section class="tp-prediction-row__arrow-icon">
        <tp-icon icon-name="arrow-right-02" weight="regular" :size="24" />
      </q-item-section>
    </q-item>
</template>

<script setup lang="ts">
import TpIcon from './TpIcon.vue'

const props = defineProps<{
  id: string
  name: string | null
  predictionType: string
  result: string
  time: string | number
  selected: boolean
}>()

const emit = defineEmits<{
  (e: 'row-click', id: string): void
  (e: 'selected-change', payload: { id: string; selected: boolean }): void
}>()

function onRowClick () {
  emit('row-click', props.id)
}
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-prediction-row {
  border-radius: 14px;
  background: var(--glass-background-light);
  transition: all 0.2s ease;
  padding: 4px 8px;
  -webkit-backdrop-filter: blur(var(--glass-blur));
  backdrop-filter: blur(var(--glass-blur));
  border: 1px solid var(--glass-border);

  &:hover {
    background: var(--glass-background);
    box-shadow: 0 4px 16px rgba(0, 0, 0, 0.06);
    border-color: var(--stroke-extra-light);
  }

  :deep(.q-focus-helper) {
    display: none;
  }

  &__prediction-type {
    max-width: 240px;
    color: var(--text-medium);
  }

  &__result {
    max-width: 160px;
  }

  &__time {
    max-width: 180px;
  }

  &__arrow-icon {
    max-width: fit-content;
  }
}
</style>
