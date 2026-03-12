<template>
  <q-page class="tp-page">
    <slot> </slot>
  </q-page>
</template>

<script setup lang="ts">

</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins' as *;

.tp-page {
  display: flex;
  flex-direction: column;
  padding: 50px 72px;
  max-width: 1440px;
  margin-inline: auto;
  width: 100%;
  min-height: 0;

  @include down(lg) {
    padding: 40px 48px;
  }

  @include down(md) {
    padding: 32px 24px;
  }

  @include down(sm) {
    padding: 24px 16px;
  }
}
</style>
