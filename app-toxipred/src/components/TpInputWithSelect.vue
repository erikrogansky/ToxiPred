<template>
  <q-input
    outlined
    v-model="inputValue"                
    @keyup.enter="$emit('enter')"
  >
    <template #append>
      <q-select
        outlined
        dense
        options-dense
        :options="options"
        option-label="label"            
        option-value="value"
        emit-value
        map-options
        style="width: 120px"
        v-model="selectedValue"
      ></q-select>
    </template>
  </q-input>
</template>

<script setup lang="ts">
const inputValue = defineModel<string>('inputValue')
const selectedValue = defineModel<string>('selectedValue')

defineProps<{
  options: Array<{ label: string; value: string }>;
}>();

defineEmits<{
  (e: 'enter'): void;
}>();
</script>

<style scoped lang="scss">
/* style away */
</style>
