<template>
  <q-layout view="hHh lpR fFf" class="tp-layout">
    <tp-header />
    <q-page-container>
      <router-view />
    </q-page-container>
    <tp-footer />
  </q-layout>
</template>

<script setup lang="ts">
import TpHeader from 'components/TpHeader.vue';
import TpFooter from 'components/TpFooter.vue';
</script>

<style scoped lang="scss">
.tp-layout {
  background: var(--page-background);
  background-attachment: fixed;
  transition: all 0.3s fade-in-out;
}
</style>
