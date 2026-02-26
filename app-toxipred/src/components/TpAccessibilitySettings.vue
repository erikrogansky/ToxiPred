<template>
  <div class="tp-accessibility-settings">
    <!-- Dark Mode Dropdown -->
    <div class="accessibility-dropdown">
      <span class="accessibility-dropdown__label">Theme</span>
      <div class="accessibility-chip">
        <span class="accessibility-chip__label">{{ darkModeModel?.label || 'Theme' }}</span>
        <tp-icon icon-name="arrow-down-02" weight="regular" :size="12" class="accessibility-chip__icon" />
        <q-menu class="accessibility-menu">
          <div class="accessibility-menu__content">
            <q-item
              v-for="option in darkModeOptions"
              :key="option.label"
              clickable
              v-close-popup
              @click="onDarkModeChange(option)"
              :active="darkModeModel?.value === option.value"
              class="accessibility-menu__item"
            >
              <q-item-section>
                <q-item-label class="paragraph-small bold">{{ option.label }}</q-item-label>
              </q-item-section>
            </q-item>
          </div>
        </q-menu>
      </div>
    </div>

    <!-- High Contrast Dropdown -->
    <div class="accessibility-dropdown">
      <span class="accessibility-dropdown__label">Contrast</span>
      <div class="accessibility-chip">
        <span class="accessibility-chip__label">{{ highContrastModel?.label || 'Contrast' }}</span>
        <tp-icon icon-name="arrow-down-02" weight="regular" :size="12" class="accessibility-chip__icon" />
        <q-menu class="accessibility-menu">
          <div class="accessibility-menu__content">
            <q-item
              v-for="option in contrastOptions"
              :key="option.label"
              clickable
              v-close-popup
              @click="onContrastChange(option)"
              :active="highContrastModel?.value === option.value"
              class="accessibility-menu__item"
            >
              <q-item-section>
                <q-item-label class="paragraph-small bold">{{ option.label }}</q-item-label>
              </q-item-section>
            </q-item>
          </div>
        </q-menu>
      </div>
    </div>

    <!-- Font Selection Dropdown -->
    <div class="accessibility-dropdown">
      <span class="accessibility-dropdown__label">Font</span>
      <div class="accessibility-chip">
        <span class="accessibility-chip__label font-preview" :class="`font-preview--${fontModel?.value || 'default'}`">{{ fontModel?.label || 'Font' }}</span>
        <tp-icon icon-name="arrow-down-02" weight="regular" :size="12" class="accessibility-chip__icon" />
        <q-menu class="accessibility-menu">
          <div class="accessibility-menu__content">
            <q-item
              v-for="option in fontOptions"
              :key="option.label"
              clickable
              v-close-popup
              @click="onFontChange(option)"
              :active="fontModel?.value === option.value"
              class="accessibility-menu__item"
            >
              <q-item-section>
                <q-item-label class="paragraph-small bold font-preview" :class="`font-preview--${option.value}`">{{ option.label }}</q-item-label>
              </q-item-section>
            </q-item>
          </div>
        </q-menu>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { storeToRefs } from 'pinia';
import { ref, watch } from 'vue';
import TpIcon from './TpIcon.vue';
import { useSettingsStore, type FontFamily } from 'src/stores/settings-store';

const settingsStore = useSettingsStore();
const { darkMode, highContrast, fontFamily } = storeToRefs(settingsStore);

// Dark mode options
const darkModeOptions = [
  { label: 'Light', value: false },
  { label: 'Dark', value: true },
];

const darkModeModel = ref(darkModeOptions.find(opt => opt.value === darkMode.value));

// Contrast options
const contrastOptions = [
  { label: 'Normal', value: false },
  { label: 'High', value: true },
];

const highContrastModel = ref(contrastOptions.find(opt => opt.value === highContrast.value));

// Font options
const fontOptions = [
  { label: 'Default', value: 'default' as FontFamily },
  { label: 'Comic Sans', value: 'comic-sans' as FontFamily },
  { label: 'OpenDyslexic', value: 'open-dyslexic' as FontFamily },
];

const fontModel = ref(fontOptions.find(opt => opt.value === fontFamily.value));

// Handlers
const onDarkModeChange = (option: { label: string; value: boolean }) => {
  darkModeModel.value = option;
  settingsStore.setDarkMode(option.value);
};

const onContrastChange = (option: { label: string; value: boolean }) => {
  highContrastModel.value = option;
  settingsStore.setHighContrast(option.value);
};

const onFontChange = (option: { label: string; value: FontFamily }) => {
  fontModel.value = option;
  settingsStore.setFontFamily(option.value);
};

// Watch for external changes
watch(darkMode, (newValue) => {
  darkModeModel.value = darkModeOptions.find(opt => opt.value === newValue);
});

watch(highContrast, (newValue) => {
  highContrastModel.value = contrastOptions.find(opt => opt.value === newValue);
});

watch(fontFamily, (newValue) => {
  fontModel.value = fontOptions.find(opt => opt.value === newValue);
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/_mixins.scss' as *;
@use 'src/css/colors/primitives' as *;

.tp-accessibility-settings {
  display: flex;
  flex-direction: column;
  gap: 8px;
  min-width: 140px;
}

.accessibility-dropdown {
  position: relative;
  display: flex;
  flex-direction: column;
  gap: 4px;

  &__label {
    font-size: 11px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--text-medium);
    padding-left: 2px;
  }
}

.accessibility-chip {
  background: var(--glass-background-light);
  border: 1px solid var(--glass-border);
  color: var(--text-brand-regular);
  cursor: pointer;
  font-size: 12px;
  font-weight: 600;
  padding: 7px 12px;
  border-radius: 20px;
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 8px;
  user-select: none;
  -webkit-backdrop-filter: blur(var(--glass-blur));
  backdrop-filter: blur(var(--glass-blur));

  &:hover {
    background: var(--glass-background);
    border-color: var(--stroke-brand-regular);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.06);
  }

  &__label {
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    flex: 1;
  }

  &__icon {
    flex-shrink: 0;
    color: var(--text-brand-regular);
  }
}

.font-preview {
  &--default {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif !important;
  }

  &--comic-sans {
    font-family: 'Comic Sans MS', 'Comic Sans', cursive !important;
  }

  &--open-dyslexic {
    font-family: 'OpenDyslexic', sans-serif !important;
  }
}
</style>

<style lang="scss">
// Global styles for the menu (not scoped)
@use 'src/css/helpers/_mixins.scss' as *;

.accessibility-menu {
  .q-menu__content {
    background: var(--glass-background) !important;
    border: 1px solid var(--glass-border) !important;
    border-radius: 16px !important;
    box-shadow: var(--glass-shadow) !important;
    -webkit-backdrop-filter: blur(var(--glass-blur-strong)) !important;
    backdrop-filter: blur(var(--glass-blur-strong)) !important;
    padding: 0 !important;
  }

  .accessibility-menu__content {
    display: flex;
    flex-direction: column;
    min-width: 140px;
  }

  .accessibility-menu__item {
    border-radius: 0;
    padding: 10px 14px;
    transition: all 0.15s ease;
    color: var(--text);

    &:first-child {
      border-radius: 16px 16px 0 0;
    }

    &:last-child {
      border-radius: 0 0 16px 16px;
    }

    &:only-child {
      border-radius: 16px;
    }

    &:hover {
      background-color: var(--surface-gray-light);
    }

    &.q-item--active {
      background-color: var(--surface-brand-extra-light);
      color: var(--text-brand-regular);
      font-weight: 700;
    }

    .q-item__label {
      color: inherit;
    }
  }
}
</style>
