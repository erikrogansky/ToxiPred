<template>
  <q-menu class="tp-menu" v-bind="$attrs">
    <div class="tp-menu__content">
      <slot />
    </div>
  </q-menu>
</template>

<script setup lang="ts">
defineOptions({ inheritAttrs: false });
</script>

<style lang="scss">
@use 'src/css/helpers/glass' as *;

.tp-menu {
  @include glass-menu;

  .tp-menu__content {
    display: flex;
    flex-direction: column;
    min-width: 160px;
  }

  .q-item {
    @include glass-menu-item;
  }
}
</style>
