<template>
  <q-page class="tp-page">
    <slot> </slot>
  </q-page>
</template>

<script setup lang="ts">

</script>

<style scoped lang="scss">
.tp-page {
  display: flex;
  flex-direction: column;
  padding: 50px 72px;
  max-width: 1440px;
  margin-inline: auto;
  width: 100%;
  min-height: 0; // Prevent unwanted stretching
}
</style>
