<template>
  <tp-page class="tp-home-page">
    <tp-button label="Click Me"></tp-button>
    <tp-button label="Click Me" variant="outline"></tp-button>
    <tp-button label="Click Me" variant="link"></tp-button>
    <h1 class="tp-main-heading">Predict dermatological toxicity from structure</h1>
    <h2 class="tp-h3 bold">OECD-aligned QSAR models with applicability domain & confidence.</h2>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpButton from 'components/TpButton.vue';
</script>

<style scoped lang="scss">
.tp-home-page {
  padding-block: 128px;

  .tp-main-heading {
    margin-bottom: 16px;
  }
}
</style>
