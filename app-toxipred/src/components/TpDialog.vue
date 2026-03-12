<template>
  <q-dialog :model-value="modelValue ?? false" @update:model-value="emit('update:modelValue', $event)" :persistent="persistent">
    <q-card class="tp-dialog">
      <q-card-section class="tp-dialog__header">
        <h3 class="tp-h3">{{ title }}</h3>
        <tp-icon-button
          icon-name="close-square"
          weight="regular"
          :size="24"
          @click="emit('close')"
        />
      </q-card-section>

      <q-card-section class="tp-dialog__content">
        <!-- Loading State -->
        <div v-if="loading" class="tp-dialog__loading">
          <q-spinner-dots size="40px" color="primary" />
          <p>{{ loadingText }}</p>
        </div>

        <!-- Error State -->
        <div v-else-if="error" class="tp-dialog__error">
          <tp-icon icon-name="warning-2" weight="regular" :size="48" />
          <p>{{ error }}</p>
          <div class="tp-dialog__actions">
            <slot name="error-actions">
              <tp-button label="Try Again" variant="outline" @click="emit('retry')" />
            </slot>
          </div>
        </div>

        <!-- Default Content -->
        <slot v-else />
      </q-card-section>
    </q-card>
  </q-dialog>
</template>

<script setup lang="ts">
import TpButton from './TpButton.vue';
import TpIconButton from './TpIconButton.vue';
import TpIcon from './TpIcon.vue';

withDefaults(defineProps<{
  modelValue?: boolean | undefined;
  title: string;
  persistent?: boolean;
  loading?: boolean;
  loadingText?: string;
  error?: string | null;
}>(), {
  persistent: false,
  loading: false,
  loadingText: 'Loading...',
  error: null,
});

const emit = defineEmits<{
  (e: 'update:modelValue', value: boolean): void;
  (e: 'close'): void;
  (e: 'retry'): void;
}>();
</script>

<style scoped lang="scss">
@use 'src/css/helpers/glass' as *;

.tp-dialog {
  @include glass-dialog;

  &__header {
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 20px 24px;
  }

  &__content {
    padding: 0 24px 24px;
  }

  &__loading {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 32px 0;
    gap: 16px;

    p {
      font-size: 14px;
      color: var(--text-medium);
    }
  }

  &__error {
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 24px 0;
    gap: 16px;
    color: var(--text-negative);

    p {
      font-size: 14px;
      text-align: center;
    }
  }

  &__actions {
    display: flex;
    justify-content: flex-end;
    gap: 12px;
    margin-top: 4px;
  }
}
</style>
