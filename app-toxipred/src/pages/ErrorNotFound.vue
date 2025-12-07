<template>
  <tp-page class="tp-404 column items-center text-center justify-center">
    <h1 class="tp-404__title">404</h1>

    <h2 class="tp-404__subtitle tp-h1 bold">Oops. Nothing here...</h2>

    <tp-button
      href="/"
      label="Go Home"
      size="big"
    />
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpButton from 'components/TpButton.vue';
</script>

<style scoped lang="scss">
.tp-404 {
  &__title {
    font-size: 80px;
  }

  &__subtitle {
    margin-bottom: 24px;
  }
}
</style>