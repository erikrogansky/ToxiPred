<template>
  <q-page class="tp-draw-page">
    <!-- ═══ Top toolbar ═══ -->
    <header class="tp-draw-toolbar">
      <div class="tp-draw-toolbar__left">
        <tp-icon-button
          icon-name="arrow-left-01"
          weight="regular"
          :size="20"
          @click="goBack"
        />
        <h3 class="tp-h4 bold">Draw Molecule</h3>
      </div>

      <div class="tp-draw-toolbar__center" v-if="currentSmiles">
        <div class="tp-draw-smiles">
          <span class="tp-draw-smiles__label">SMILES</span>
          <code class="tp-draw-smiles__value">{{ currentSmiles }}</code>
          <tp-icon-button
            icon-name="copy"
            weight="regular"
            :size="16"
            @click="copySmiles"
            title="Copy SMILES"
          />
        </div>
      </div>

      <div class="tp-draw-toolbar__right">
        <div class="tp-draw-model-select">
          <label class="tp-draw-model-select__label">Model</label>
          <div class="tp-draw-model-select__chip" @click.stop>
            <span>{{ currentModelLabel }}</span>
            <tp-icon icon-name="arrow-down-02" weight="regular" :size="14" />
            <q-menu class="tp-draw-model-menu">
              <div class="tp-draw-model-menu__content">
                <q-item
                  v-for="opt in modelOptions"
                  :key="opt.value"
                  clickable
                  v-close-popup
                  @click="selectedModel = opt.value"
                  :active="selectedModel === opt.value"
                  class="tp-draw-model-menu__item"
                >
                  <q-item-section>
                    <q-item-label class="paragraph-small bold">{{ opt.label }}</q-item-label>
                  </q-item-section>
                </q-item>
              </div>
            </q-menu>
          </div>
        </div>

        <tp-button
          label="Predict"
          :disabled="!currentSmiles || !selectedModel || predicting"
          @click="submitPrediction"
        />
      </div>
    </header>

    <!-- ═══ Ketcher editor ═══ -->
    <div class="tp-draw-editor">
      <tp-ketcher-editor
        ref="editorRef"
        :initial-smiles="initialSmiles"
        @ready="onKetcherReady"
      />
    </div>

    <!-- ═══ Snackbar for copy feedback ═══ -->
    <Transition name="fade">
      <div v-if="showCopied" class="tp-draw-toast">
        SMILES copied to clipboard
      </div>
    </Transition>
  </q-page>
</template>

<script setup lang="ts">
import { ref, computed, onBeforeUnmount } from 'vue';
import { useRouter, useRoute } from 'vue-router';
import TpKetcherEditor from 'components/TpKetcherEditor.vue';
import TpButton from 'components/TpButton.vue';
import TpIcon from 'components/TpIcon.vue';
import TpIconButton from 'components/TpIconButton.vue';
import { useModelsStore, type TestType, type PredictionTarget } from 'src/stores/models-store';
import { useJobsStore } from 'src/stores/jobs-store';
import { api } from 'src/boot/axios';

const router = useRouter();
const route = useRoute();
const modelsStore = useModelsStore();
const jobsStore = useJobsStore();

// ─── State ─────────────────────────────────────────
const editorRef = ref<InstanceType<typeof TpKetcherEditor>>();
const currentSmiles = ref('');
const selectedModel = ref((route.query.model as string) || '');
const initialSmiles = (route.query.smiles as string) || '';
const predicting = ref(false);
const showCopied = ref(false);

let smilesPollerHandle: ReturnType<typeof setInterval> | null = null;

// ─── Model options ─────────────────────────────────
const testTypeLabels: Record<TestType, string> = {
  in_vitro: 'In Vitro',
  in_vivo: 'In Vivo',
  in_chemico: 'In Chemico',
};
const predictionTargetLabels: Record<PredictionTarget, string> = {
  photo_irritation: 'Photo Irritation',
  photo_toxicity: 'Phototoxicity',
};

const modelOptions = computed(() =>
  modelsStore.getModels.map((name) => {
    const d = modelsStore.getModelDetail(name);
    const parts: string[] = [];
    if (d?.prediction_target) parts.push(predictionTargetLabels[d.prediction_target] || d.prediction_target);
    if (d?.test_type) parts.push(testTypeLabels[d.test_type] || d.test_type);
    return { label: parts.length ? parts.join(' – ') : name, value: name };
  }),
);

const currentModelLabel = computed(() => {
  const found = modelOptions.value.find((o) => o.value === selectedModel.value);
  return found?.label || 'Select model';
});

// ─── Ketcher callbacks ─────────────────────────────
// eslint-disable-next-line @typescript-eslint/no-explicit-any
function onKetcherReady(ketcher: any) {
  // Poll SMILES every 1.5s (Ketcher doesn't expose reliable change events in all versions)
  smilesPollerHandle = setInterval(() => {
    void (async () => {
      try {
        const smiles = await ketcher.getSmiles() as string;
        // Ketcher returns empty or single-element SMILES for blank canvas
        currentSmiles.value = smiles && smiles.length > 1 ? smiles : '';
      } catch {
        /* ignore transient errors during edits */
      }
    })();
  }, 1500);
}

// ─── Actions ───────────────────────────────────────
function goBack() {
  if (window.history.length > 1) {
    router.back();
  } else {
    void router.push({ name: 'home' });
  }
}

async function copySmiles() {
  if (!currentSmiles.value) return;
  await navigator.clipboard.writeText(currentSmiles.value);
  showCopied.value = true;
  setTimeout(() => (showCopied.value = false), 2000);
}

async function submitPrediction() {
  if (!currentSmiles.value || !selectedModel.value) return;
  predicting.value = true;

  try {
    const response = await api.post(`/jobs/predict/${selectedModel.value}`, {
      query: currentSmiles.value,
    });

    if (response.data.error) {
      alert(`Error: ${response.data.error}`);
      return;
    }
    if (!response.data.job_id) {
      alert('Error: No job ID returned from server');
      return;
    }

    jobsStore.upsert({
      id: response.data.job_id,
      model: selectedModel.value,
      state: 'PENDING',
      percent: 0,
    });
    jobsStore.startTracking(response.data.job_id);

    await router.push({ name: 'workspace' });
  } catch {
    alert('Error submitting prediction. Please try again.');
  } finally {
    predicting.value = false;
  }
}

// ─── Cleanup ───────────────────────────────────────
onBeforeUnmount(() => {
  if (smilesPollerHandle) {
    clearInterval(smilesPollerHandle);
    smilesPollerHandle = null;
  }
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

/* ═══ Page layout ═══ */
.tp-draw-page {
  display: flex;
  flex-direction: column;
  padding: 0 !important;
  max-width: 1440px;
  margin-inline: auto;
  width: 100%;
  /* Override Quasar's computed min-height so the page doesn't push beyond the viewport */
  min-height: 0 !important;
  height: calc(100dvh - 60px);
  max-height: calc(100dvh - 60px);
  overflow: hidden;
}

/* ═══ Toolbar ═══ */
.tp-draw-toolbar {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
  padding: 12px 32px;
  border-bottom: 1px solid var(--glass-border);
  background: var(--glass-background);
  backdrop-filter: blur(var(--glass-blur));
  -webkit-backdrop-filter: blur(var(--glass-blur));
  flex-shrink: 0;
  z-index: 10;

  @include down(md) {
    flex-wrap: wrap;
    padding: 10px 16px;
  }

  &__left {
    display: flex;
    align-items: center;
    gap: 12px;
    flex-shrink: 0;

    h3 {
      white-space: nowrap;
    }
  }

  &__center {
    flex: 1;
    min-width: 0;
    display: flex;
    justify-content: center;

    @include down(md) {
      order: 3;
      flex-basis: 100%;
      justify-content: flex-start;
    }
  }

  &__right {
    display: flex;
    align-items: center;
    gap: 12px;
    flex-shrink: 0;
  }
}

/* ═══ SMILES display chip ═══ */
.tp-draw-smiles {
  display: flex;
  align-items: center;
  gap: 8px;
  background: var(--glass-background-light);
  border: 1px solid var(--glass-border);
  border-radius: 28px;
  padding: 6px 12px 6px 16px;
  max-width: 480px;
  overflow: hidden;

  &__label {
    font-size: 11px;
    font-weight: 900;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    color: var(--text-brand-regular);
    flex-shrink: 0;
  }

  &__value {
    font-size: 13px;
    color: var(--text);
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    font-family: 'SF Mono', 'Fira Code', monospace;
  }
}

/* ═══ Model selector ═══ */
.tp-draw-model-select {
  display: flex;
  align-items: center;
  gap: 8px;

  &__label {
    font-size: 13px;
    font-weight: 600;
    color: var(--text-medium);
    white-space: nowrap;

    @include down(sm) {
      display: none;
    }
  }

  &__chip {
    background: var(--surface-brand-extra-light);
    border: 1px solid var(--stroke-brand-regular);
    color: var(--text-brand-regular);
    cursor: pointer;
    font-size: 13px;
    font-weight: 600;
    padding: 6px 12px;
    border-radius: 18px;
    display: flex;
    align-items: center;
    gap: 6px;
    user-select: none;
    transition: all 0.2s ease;

    &:hover {
      background: var(--surface-brand-light);
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.06);
    }

    span {
      max-width: 200px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
  }
}

.tp-draw-model-menu {
  background: var(--glass-background) !important;
  border: 1px solid var(--glass-border) !important;
  border-radius: 16px !important;
  box-shadow: var(--glass-shadow) !important;
  backdrop-filter: blur(var(--glass-blur-strong)) !important;
  -webkit-backdrop-filter: blur(var(--glass-blur-strong)) !important;

  &__content {
    display: flex;
    flex-direction: column;
    min-width: 220px;
  }

  &__item {
    border-radius: 10px;
    padding: 10px 12px;
    transition: all 0.15s ease;

    &:hover {
      background: var(--surface-brand-extra-light);
    }

    &.q-item--active {
      background: var(--surface-brand-extra-light);
      color: var(--text-brand-regular);
      font-weight: 600;
    }
  }
}

/* ═══ Editor container ═══ */
.tp-draw-editor {
  flex: 1;
  min-height: 0;
  position: relative;
  overflow: hidden;

  :deep(.tp-ketcher-editor) {
    height: 100%;
    overflow: hidden;
  }
}

/* ═══ Toast notification ═══ */
.tp-draw-toast {
  position: fixed;
  bottom: 24px;
  left: 50%;
  transform: translateX(-50%);
  background: var(--surface-brand-regular);
  color: var(--text-inverse);
  font-size: 13px;
  font-weight: 600;
  padding: 10px 20px;
  border-radius: 360px;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
  z-index: 9999;
  pointer-events: none;
}

/* ═══ Transitions ═══ */
.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.25s ease;
}
.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}
</style>
