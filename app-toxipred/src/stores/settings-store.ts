import { defineStore, acceptHMRUpdate } from 'pinia';

export const useSettingsStore = defineStore('settings', {
  state: () => ({
    darkMode: true,
  }),

  getters: {
    //getDarkMode: (state) => state.darkMode,
  },

  actions: {
    setDarkMode(darkMode: boolean) {
        this.darkMode = darkMode
    },

    toggleDarkMode() {
        this.darkMode = !this.darkMode
    }
  },
});

if (import.meta.hot) {
  import.meta.hot.accept(acceptHMRUpdate(useSettingsStore, import.meta.hot));
}
