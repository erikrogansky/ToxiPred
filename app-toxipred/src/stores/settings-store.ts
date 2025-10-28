import { defineStore, acceptHMRUpdate } from 'pinia';

export const useSettingsStore = defineStore('settings', {
  state: () => ({
    colorMode: 'dark'
  }),

  getters: {
    colorMode: (state) => state.colorMode,
  },

  actions: {
    setColorMode(mode: string) {
        this.colorMode = mode
    }
  },
});

if (import.meta.hot) {
  import.meta.hot.accept(acceptHMRUpdate(useSettingsStore, import.meta.hot));
}
