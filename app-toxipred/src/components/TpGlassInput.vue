<template>
  <div v-if="type === 'password'" class="tp-glass-input-wrapper tp-glass-input-wrapper--password">
    <input
      :value="modelValue"
      @input="$emit('update:modelValue', ($event.target as HTMLInputElement).value)"
      :type="passwordVisible ? 'text' : 'password'"
      :placeholder="placeholder"
      :readonly="readonly"
      class="tp-glass-input-wrapper__field"
      @keyup.enter="$emit('enter')"
    />
    <button class="tp-glass-input-wrapper__toggle" @click="passwordVisible = !passwordVisible" type="button">
      <tp-icon :icon-name="passwordVisible ? 'eye-slash' : 'eye'" weight="regular" :size="18" />
    </button>
  </div>

  <input
    v-else
    :value="modelValue"
    @input="$emit('update:modelValue', ($event.target as HTMLInputElement).value)"
    :type="type"
    :placeholder="placeholder"
    :readonly="readonly"
    class="tp-glass-input"
    :class="{ 'tp-glass-input--monospace': monospace }"
    @keyup.enter="$emit('enter')"
  />
</template>

<script setup lang="ts">
import { ref } from 'vue';
import TpIcon from './TpIcon.vue';

withDefaults(defineProps<{
  modelValue?: string;
  placeholder?: string;
  type?: 'text' | 'password';
  readonly?: boolean;
  monospace?: boolean;
}>(), {
  modelValue: '',
  placeholder: '',
  type: 'text',
  readonly: false,
  monospace: false,
});

defineEmits<{
  (e: 'update:modelValue', value: string): void;
  (e: 'enter'): void;
}>();

const passwordVisible = ref(false);
</script>

<style scoped lang="scss">
@use 'src/css/helpers/glass' as *;

.tp-glass-input {
  @include glass-input;

  &--monospace {
    font-family: monospace;
    letter-spacing: 0.5px;
  }
}

.tp-glass-input-wrapper--password {
  @include glass-password-wrapper;
}

.tp-glass-input-wrapper__field {
  flex: 1;
  border: none;
  outline: none;
  background: transparent;
  backdrop-filter: none;
  padding: 12px 16px;
  font-size: 14px;
  color: var(--text);

  &::placeholder {
    color: var(--text-medium);
  }
}

.tp-glass-input-wrapper__toggle {
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 6px;
  border: none;
  background: transparent;
  color: var(--text-light);
  cursor: pointer;
  border-radius: 6px;
  transition: all 0.15s ease;

  &:hover {
    color: var(--text);
  }
}
</style>
