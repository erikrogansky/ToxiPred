<template>
  <div ref="ketcherHost" class="tp-ketcher-editor">
    <div v-if="loading" class="tp-ketcher-loading">
      <q-spinner-dots size="48px" color="primary" />
      <p class="paragraph-small">Loading molecular editor…</p>
    </div>
    <div v-if="error" class="tp-ketcher-error">
      <p class="paragraph-small">{{ error }}</p>
      <tp-button label="Retry" variant="outline" size="small" @click="initKetcher" />
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, onBeforeUnmount } from 'vue';
import TpButton from './TpButton.vue';

const ketcherHost = ref<HTMLDivElement>();
const loading = ref(true);
const error = ref<string | null>(null);

// eslint-disable-next-line @typescript-eslint/no-explicit-any
let reactRoot: any = null;
// eslint-disable-next-line @typescript-eslint/no-explicit-any
let ketcherInstance: any = null;
let ketcherReady = false;

const emit = defineEmits<{
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  (e: 'ready', ketcher: any): void;
}>();

const props = defineProps<{
  initialSmiles?: string;
}>();

async function initKetcher() {
  if (!ketcherHost.value) return;

  loading.value = true;
  error.value = null;

  try {
    const [React, ReactDOM, ketcherReact] = await Promise.all([
      import('react'),
      import('react-dom/client'),
      import('ketcher-react'),
    ]);

    // @ts-expect-error ketcher-standalone package.json exports missing types condition
    const ketcherStandalone = await import('ketcher-standalone');

    // Load Ketcher CSS (Vite handles this as a side-effect import — injected once globally)
    await import('ketcher-react/dist/index.css');

    if (!ketcherHost.value) return;

    const structServiceProvider = new ketcherStandalone.StandaloneStructServiceProvider();

    reactRoot = ReactDOM.createRoot(ketcherHost.value);
    reactRoot.render(
      React.createElement(ketcherReact.Editor, {
        staticResourcesUrl: '',
        structServiceProvider,
        errorHandler: (msg: string) => console.warn('[Ketcher]', msg),
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        onInit: (ketcher: any) => {
          ketcherInstance = ketcher;
          ketcherReady = true;
          loading.value = false;

          if (props.initialSmiles) {
            ketcher.setMolecule(props.initialSmiles).catch(() => {/* ignore */});
          }

          emit('ready', ketcher);
        },
      })
    );
  } catch (e) {
    console.error('Failed to load Ketcher:', e);
    error.value = 'Failed to load the molecular editor. Please try again.';
    loading.value = false;
  }
}

onMounted(() => {
  void initKetcher();
});

onBeforeUnmount(() => {
  if (reactRoot) {
    reactRoot.unmount();
    reactRoot = null;
  }
  ketcherInstance = null;
  ketcherReady = false;
});

async function getSmiles(): Promise<string> {
  if (!ketcherInstance) return '';
  try {
    return await ketcherInstance.getSmiles();
  } catch {
    return '';
  }
}

async function getMolfile(): Promise<string> {
  if (!ketcherInstance) return '';
  try {
    return await ketcherInstance.getMolfile();
  } catch {
    return '';
  }
}

async function setMolecule(mol: string): Promise<void> {
  if (!ketcherInstance) return;
  try {
    await ketcherInstance.setMolecule(mol);
  } catch (e) {
    console.error('Failed to set molecule:', e);
  }
}

defineExpose({
  getSmiles,
  getMolfile,
  setMolecule,
  getKetcher: () => ketcherInstance,
  isReady: () => ketcherReady,
});
</script>

<style scoped lang="scss">
.tp-ketcher-editor {
  width: 100%;
  height: 100%;
  position: relative;
  min-height: 400px;
}

.tp-ketcher-loading,
.tp-ketcher-error {
  position: absolute;
  inset: 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: 16px;
  z-index: 10;
  background: var(--glass-background);
}
</style>
