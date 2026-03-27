<template>
  <tp-page class="tp-demos-page column items-center">

    <!-- ========== HERO ========== -->
    <section class="tp-demos-hero column items-center">
      <h1 class="tp-main-heading">Demo Predictions</h1>
      <p class="paragraph tp-demos-hero__subtitle">
        Explore pre-computed predictions to see ToxiPred in action. Each demo showcases a real compound
        with full results — including feature importance, atom contributions, and applicability domain assessment.
      </p>
    </section>

    <!-- ========== LOADING ========== -->
    <div v-if="loading" class="column items-center q-pa-xl">
      <q-spinner-dots size="48px" color="primary" />
      <p class="q-mt-md paragraph">Loading demos...</p>
    </div>

    <!-- ========== ERROR ========== -->
    <div v-else-if="error" class="column items-center q-pa-xl">
      <tp-icon icon-name="warning-2" weight="regular" :size="48" />
      <p class="q-mt-md paragraph">{{ error }}</p>
      <tp-button label="Retry" variant="outline" @click="fetchDemos" class="q-mt-md" />
    </div>

    <!-- ========== EMPTY ========== -->
    <div v-else-if="demos.length === 0" class="column items-center q-pa-xl">
      <tp-icon icon-name="discover" weight="regular" :size="48" />
      <p class="q-mt-md paragraph">No demo predictions available yet.</p>
    </div>

    <!-- ========== DEMO CARDS ========== -->
    <section v-else class="tp-demos-grid">
      <div v-for="demo in demos" :key="demo.job_id" class="tp-demo-card">
        <div class="tp-demo-card__badge">{{ demo.model }}</div>
        <h3 class="tp-demo-card__title">{{ demo.demo_title || demo.trivial_name || demo.name || 'Demo' }}</h3>
        <p class="paragraph-small tp-demo-card__description">
          {{ demo.demo_description || 'Explore this prediction to see the full results, feature importance, and atom contributions.' }}
        </p>

        <div class="tp-demo-card__details">
          <div class="tp-demo-card__detail" v-if="demo.trivial_name">
            <span class="tp-demo-card__detail-label">Compound</span>
            <span class="tp-demo-card__detail-value">{{ demo.trivial_name }}</span>
          </div>
          <div class="tp-demo-card__detail" v-if="demo.formula">
            <span class="tp-demo-card__detail-label">Formula</span>
            <span class="tp-demo-card__detail-value">{{ demo.formula }}</span>
          </div>
          <div class="tp-demo-card__detail" v-if="demo.prediction != null">
            <span class="tp-demo-card__detail-label">Result</span>
            <span class="tp-demo-card__detail-value tp-demo-card__detail-value--result" :class="predictionClass(demo)">
              {{ predictionLabel(demo) }}
            </span>
          </div>
        </div>

        <tp-button label="View Prediction" @click="openDemo(demo.job_id)" />
      </div>
    </section>

    <!-- ========== CTA ========== -->
    <section v-if="demos.length > 0" class="tp-demos-cta">
      <p class="paragraph">Ready to try your own compound?</p>
      <tp-button label="Start predicting" @click="$router.push('/')" />
    </section>

  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpButton from 'components/TpButton.vue';
import TpIcon from 'components/TpIcon.vue';
import { ref, onMounted } from 'vue';
import { useRouter } from 'vue-router';
import { useModelsStore } from 'src/stores/models-store';
import { api } from 'src/boot/axios';

interface Demo {
  job_id: string;
  model: string;
  name: string | null;
  trivial_name: string | null;
  formula: string | null;
  canonical_smiles: string | null;
  demo_title: string | null;
  demo_description: string | null;
  prediction: number | number[] | null;
}

const router = useRouter();
const modelsStore = useModelsStore();

const loading = ref(true);
const error = ref<string | null>(null);
const demos = ref<Demo[]>([]);

function predictionValue(demo: Demo): number | null {
  if (demo.prediction == null) return null;
  return Array.isArray(demo.prediction) ? (demo.prediction[0] ?? null) : demo.prediction;
}

function predictionLabel(demo: Demo): string {
  const val = predictionValue(demo);
  if (val == null) return 'N/A';
  const detail = modelsStore.getModelDetail(demo.model);
  return val === 1
    ? (detail?.positive_label || 'Toxic')
    : (detail?.negative_label || 'Non-toxic');
}

function predictionClass(demo: Demo): string {
  const val = predictionValue(demo);
  return val === 1 ? 'tp-demo-card__detail-value--positive' : 'tp-demo-card__detail-value--negative';
}

function openDemo(jobId: string) {
  void router.push({ name: 'job-overview', params: { job_id: jobId } });
}

async function fetchDemos() {
  loading.value = true;
  error.value = null;
  try {
    const response = await api.get('/demos');
    demos.value = response.data;
  } catch {
    error.value = 'Failed to load demo predictions.';
  } finally {
    loading.value = false;
  }
}

onMounted(() => void fetchDemos());
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-demos-page {
  padding-block: 0;
  gap: 64px;
}

// ─── HERO ──────────────────────────────────────────
.tp-demos-hero {
  padding-top: 96px;
  text-align: center;

  .tp-main-heading {
    margin-bottom: 16px;
  }

  &__subtitle {
    max-width: 640px;
    color: var(--text-medium);
    line-height: 170%;
  }
}

// ─── DEMO CARDS ────────────────────────────────────
.tp-demos-grid {
  width: 100%;
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 28px;

  @include down(lg) {
    grid-template-columns: repeat(2, 1fr);
  }

  @include down(sm) {
    grid-template-columns: 1fr;
  }
}

.tp-demo-card {
  background: var(--glass-background);
  border: 1px solid var(--glass-border);
  border-radius: 20px;
  padding: 36px 32px;
  backdrop-filter: blur(var(--glass-blur));
  -webkit-backdrop-filter: blur(var(--glass-blur));
  display: flex;
  flex-direction: column;
  gap: 16px;
  transition: transform 0.2s ease, box-shadow 0.2s ease;

  &:hover {
    transform: translateY(-4px);
    box-shadow: var(--glass-shadow-elevated);
  }

  &__badge {
    display: inline-block;
    width: fit-content;
    font-size: 11px;
    font-weight: 900;
    text-transform: uppercase;
    letter-spacing: 0.6px;
    padding: 4px 12px;
    border-radius: 360px;
    background: var(--surface-brand-extra-light);
    color: var(--text-brand-regular);
  }

  &__title {
    margin: 0;
  }

  &__description {
    color: var(--text-medium);
    line-height: 160%;
    flex: 1;
  }

  &__details {
    display: flex;
    flex-direction: column;
    gap: 8px;
    padding: 16px 0;
    border-top: 1px solid var(--glass-border);
  }

  &__detail {
    display: flex;
    justify-content: space-between;
    align-items: center;

    &-label {
      font-size: 13px;
      color: var(--text-light);
      font-weight: 600;
    }

    &-value {
      font-size: 13px;
      font-weight: 700;

      &--result {
        padding: 2px 10px;
        border-radius: 360px;
        font-size: 12px;
      }

      &--positive {
        background: color-mix(in srgb, var(--status-error) 15%, transparent);
        color: var(--status-error);
      }

      &--negative {
        background: color-mix(in srgb, var(--status-success) 15%, transparent);
        color: var(--status-success);
      }
    }
  }
}

// ─── CTA ───────────────────────────────────────────
.tp-demos-cta {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 16px;
  padding-bottom: 96px;

  p {
    color: var(--text-medium);
  }
}
</style>
