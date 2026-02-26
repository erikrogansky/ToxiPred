import { defineStore, acceptHMRUpdate } from 'pinia';

export type FontFamily = 'default' | 'comic-sans' | 'open-dyslexic';

export const useSettingsStore = defineStore('settings', {
  state: () => ({
    darkMode: true,
    highContrast: false,
    fontFamily: 'default' as FontFamily,
  }),

  getters: {
    theme: (state): string => {
      if (state.highContrast) {
        return state.darkMode ? 'hc-dark' : 'hc-light';
      }
      return state.darkMode ? 'dark' : 'light';
    },
  },

  actions: {
    setDarkMode(darkMode: boolean) {
        this.darkMode = darkMode
    },

    toggleDarkMode() {
        this.darkMode = !this.darkMode
    },

    setHighContrast(highContrast: boolean) {
        this.highContrast = highContrast
    },

    toggleHighContrast() {
        this.highContrast = !this.highContrast
    },

    setFontFamily(fontFamily: FontFamily) {
        this.fontFamily = fontFamily
    },
  },
});

if (import.meta.hot) {
  import.meta.hot.accept(acceptHMRUpdate(useSettingsStore, import.meta.hot));
}
