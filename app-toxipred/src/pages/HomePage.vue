<template>
  <tp-page class="tp-home-page column items-center">

    <!-- ========== HERO ========== -->
    <section class="tp-hero column items-center">
      <h1 class="tp-main-heading">Predict dermatological toxicity from structure</h1>
      <h2 class="tp-h3 bold">OECD-aligned QSAR models with applicability domain &amp; confidence.</h2>

      <div class="row items-center justify-center full-width q-mt-md q-gutter-md">
        <tp-button label="See demo predictions" />
        <tp-button label="Documentation" variant="outline" href="/documentation" />
      </div>

      <div class="row items-top justify-center full-width tp-input-section">
        <tp-input-with-select select-label="Prediction" label="Enter SMILES, CAS number, or trivial name" hint="For batch submissions, separate the identifiers by a comma" v-model:input-value="inputValue" v-model:selected-value="selectedValue" :options="modelOptions" @enter="submitPrediction"/>
        <tp-button :disabled="!inputValue || !selectedValue" style="margin: 0 0 20px 16px;" @click="submitPrediction" label="Predict" variant="outline" size="regular"></tp-button>
      </div>
    </section>

    <!-- ========== FEATURE CARDS ========== -->
    <section class="tp-features">
      <h2 class="tp-section-title">Why ToxiPred?</h2>
      <div class="tp-features__grid">
        <div v-for="(card, i) in featureCards" :key="i" class="tp-feature-card">
          <div class="tp-feature-card__icon">
            <tp-icon :icon-name="card.icon" weight="regular" :size="28" />
          </div>
          <h4>{{ card.title }}</h4>
          <p class="paragraph-small">{{ card.description }}</p>
        </div>
      </div>
    </section>

    <!-- ========== HOW IT WORKS ========== -->
    <section class="tp-how-it-works">
      <h2 class="tp-section-title">How it works</h2>
      <div class="tp-steps">
        <div v-for="(step, i) in steps" :key="i" class="tp-step">
          <div class="tp-step__number">{{ i + 1 }}</div>
          <div class="tp-step__content">
            <h4>{{ step.title }}</h4>
            <p class="paragraph-small">{{ step.text }}</p>
          </div>
        </div>
      </div>
    </section>

    <!-- ========== STATS BANNER ========== -->
    <section class="tp-stats-banner">
      <div v-for="(stat, i) in stats" :key="i" class="tp-stat">
        <span class="tp-stat__number">{{ stat.value }}</span>
        <span class="tp-stat__label paragraph-small">{{ stat.label }}</span>
      </div>
    </section>

    <!-- ========== MODELS SHOWCASE ========== -->
    <section class="tp-models">
      <h2 class="tp-section-title">Available models</h2>
      <div class="tp-models__grid">
        <div v-for="(model, i) in models" :key="i" class="tp-model-card">
          <div class="tp-model-card__header">
            <span class="tp-model-card__badge">{{ model.badge }}</span>
            <h4>{{ model.name }}</h4>
          </div>
          <ul class="tp-model-card__details">
            <li v-for="(detail, j) in model.details" :key="j">
              <tp-icon icon-name="check" weight="regular" :size="16" />
              <span class="paragraph-small">{{ detail }}</span>
            </li>
          </ul>
          <tp-button :label="`Try ${model.name}`" variant="outline" size="small" @click="selectAndScroll(model.value)" />
        </div>
      </div>
    </section>

    <!-- ========== CTA BANNER ========== -->
    <section class="tp-cta-banner">
      <div class="tp-cta-banner__content">
        <h2>Ready to predict?</h2>
        <p class="paragraph">Enter a SMILES string, CAS number, or compound name above to get started — no registration required.</p>
        <div class="row q-gutter-md">
          <tp-button label="Start predicting" @click="scrollToInput" />
          <tp-button label="View documentation" variant="outline" href="/documentation" />
        </div>
      </div>
    </section>

    <!-- ========== FAQ ========== -->
    <section class="tp-faq">
      <h2 class="tp-section-title">Frequently asked questions</h2>
      <div class="tp-faq__list">
        <div v-for="(item, i) in faqItems" :key="i" class="tp-faq__item" :class="{ 'tp-faq__item--open': faqOpen === i }">
          <button class="tp-faq__question" @click="faqOpen = faqOpen === i ? -1 : i">
            <span class="bold">{{ item.q }}</span>
            <tp-icon :icon-name="faqOpen === i ? 'minus' : 'add'" weight="regular" :size="18" />
          </button>
          <Transition name="faq-expand">
            <div v-if="faqOpen === i" class="tp-faq__answer">
              <p class="paragraph-small" v-html="item.a"></p>
            </div>
          </Transition>
        </div>
      </div>
    </section>

  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpInputWithSelect from 'components/TpInputWithSelect.vue';
import TpButton from 'components/TpButton.vue';
import TpIcon from 'components/TpIcon.vue';
import { ref, computed } from 'vue';

import { useModelsStore, type TestType, type PredictionTarget } from 'src/stores/models-store';
import { useJobsStore } from 'src/stores/jobs-store';

const modelsStore = useModelsStore()

// Labels for display
const testTypeLabels: Record<TestType, string> = {
  'in_vitro': 'In Vitro',
  'in_vivo': 'In Vivo',
  'in_chemico': 'In Chemico',
}

const predictionTargetLabels: Record<PredictionTarget, string> = {
  'photo_irritation': 'Photo Irritation',
  'photo_toxicity': 'Phototoxicity',
}

// Create model options with formatted labels (Test Type + Prediction Target)
const modelOptions = computed(() => {
  return modelsStore.getModels.map(modelName => {
    const detail = modelsStore.getModelDetail(modelName)
    const parts: string[] = []
    
    if (detail?.prediction_target) {
      parts.push(predictionTargetLabels[detail.prediction_target] || detail.prediction_target)
    }
    if (detail?.test_type) {
      parts.push(testTypeLabels[detail.test_type] || detail.test_type)
    }
    
    const label = parts.length > 0 ? parts.join(' - ') : modelName
    
    return {
      label,
      value: modelName
    }
  })
})

const inputValue = ref("");
const selectedValue = ref("");

import { api } from 'src/boot/axios';
import { useRouter } from 'vue-router';
const router = useRouter();

const submitPrediction = async () => {
  try {
    const response = await api.post(`/jobs/predict/${selectedValue.value}`, {
      query: inputValue.value
    });

    if (response.data.error) {
      console.error('Prediction error:', response.data.error);
      alert(`Error: ${response.data.error}`);
      return;
    }

    if (!response.data.job_id) {
      console.error('No job_id in response:', response.data);
      alert('Error: No job ID returned from server');
      return;
    }

    const jobs = useJobsStore()
    jobs.upsert({
      id: response.data.job_id,
      model: selectedValue.value,
      state: 'PENDING',
      percent: 0,
    })
    jobs.startTracking(response.data.job_id)

    await router.push({ name: 'workspace' });

  } catch (error) {
    console.error('Error making prediction:', error);
    alert('Error submitting prediction. Please try again.');
  }
};

// ─── SECTION DATA ────────────────────────────────────────

const featureCards = [
  {
    icon: 'shield' as const,
    title: 'OECD-aligned',
    description: 'Models built following OECD principles for the validation of QSAR models, ensuring regulatory acceptance.',
  },
  {
    icon: 'chart' as const,
    title: 'Applicability domain',
    description: 'Each prediction is checked against the model\'s applicability domain so you know when results are reliable.',
  },
  {
    icon: 'discover' as const,
    title: 'Explainability',
    description: 'Feature importance and SHAP-based explanations help you understand what drives each prediction.',
  },
  {
    icon: 'task' as const,
    title: 'Batch processing',
    description: 'Submit multiple compounds at once — via comma-separated input or bulk upload — and track all jobs in your workspace.',
  },
  {
    icon: 'share' as const,
    title: 'Shareable results',
    description: 'Generate a unique link for any prediction job to share results with colleagues or include in reports.',
  },
  {
    icon: 'code' as const,
    title: 'Open & transparent',
    description: 'Descriptors, model type, and training details are fully documented. No black boxes.',
  },
]

const steps = [
  {
    title: 'Choose a model',
    text: 'Select from available QSAR models — each targeting a specific dermatological endpoint such as phototoxicity or photo-irritation.',
  },
  {
    title: 'Enter your compound',
    text: 'Provide a SMILES string, CAS registry number, or trivial name. For batch predictions, separate identifiers with commas.',
  },
  {
    title: 'Get your prediction',
    text: 'Receive a binary classification with confidence score, applicability domain assessment, and molecular descriptor breakdown.',
  },
  {
    title: 'Explore & share',
    text: 'Dive into feature importance charts, XSMILES visualisations, and share results via a unique link.',
  },
]

const stats = [
  { value: '2', label: 'QSAR models' },
  { value: '32', label: 'Molecular descriptors' },
  { value: '100%', label: 'Free to use' },
  { value: '<5s', label: 'Prediction time' },
]

const models = [
  {
    name: 'XGBoost',
    value: 'XGBoost',
    badge: 'In Vitro',
    details: [
      'Phototoxicity prediction',
      '11 molecular descriptors',
      'Tree-based SHAP explanations',
      'Applicability domain check',
    ],
  },
  {
    name: 'Ensemble',
    value: 'Ensamble',
    badge: 'In Chemico',
    details: [
      'Photo-irritation prediction',
      '21 molecular descriptors',
      'Ensemble voting classifier',
      'Applicability domain check',
    ],
  },
]

const faqItems = [
  {
    q: 'What input formats are supported?',
    a: 'You can enter <strong>SMILES strings</strong>, <strong>CAS registry numbers</strong>, or <strong>trivial (common) compound names</strong>. ToxiPred will resolve them automatically. For batch jobs, separate multiple identifiers with commas.',
  },
  {
    q: 'What does "applicability domain" mean?',
    a: 'The applicability domain defines the chemical space in which a QSAR model can make reliable predictions. If your compound falls outside this domain, ToxiPred will flag the prediction as potentially unreliable — giving you an honest confidence indicator.',
  },
  {
    q: 'Are the models validated?',
    a: 'Yes. Both models are built following <strong>OECD principles</strong> for QSAR validation, including a defined endpoint, an unambiguous algorithm, a defined applicability domain, appropriate measures of goodness-of-fit, and a mechanistic interpretation where possible.',
  },
  {
    q: 'Do I need to create an account?',
    a: 'No. ToxiPred is <strong>completely free</strong> and requires no registration. Simply enter your compound and get results instantly.',
  },
  {
    q: 'Can I share my results?',
    a: 'Yes. Every prediction job gets a <strong>unique shareable link</strong> that you can send to colleagues or reference in publications and regulatory submissions.',
  },
  {
    q: 'What dermatological endpoints are covered?',
    a: 'Currently ToxiPred covers <strong>phototoxicity</strong> (in vitro, XGBoost model) and <strong>photo-irritation</strong> (in chemico, ensemble model). Additional endpoints may be added in the future.',
  },
]

const faqOpen = ref(-1)

function scrollToInput() {
  const inputEl = document.querySelector('.tp-input-section')
  inputEl?.scrollIntoView({ behavior: 'smooth', block: 'center' })
}

function selectAndScroll(modelValue: string) {
  selectedValue.value = modelValue
  scrollToInput()
}
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

// ─── SHARED ────────────────────────────────────────
.tp-section-title {
  text-align: center;
  margin-bottom: 48px;
}

.tp-home-page {
  padding-block: 0;
  gap: 96px;

  @include down(md) {
    gap: 64px;
  }
}

// ─── HERO ──────────────────────────────────────────
.tp-hero {
  padding-block: 128px 0;
  width: 100%;

  .tp-main-heading {
    margin-bottom: 16px;
    text-align: center;
  }

  .tp-h3 {
    text-align: center;
  }

  .tp-input-section {
    margin-top: 72px;
  }
}

// ─── FEATURE CARDS ─────────────────────────────────
.tp-features {
  width: 100%;

  &__grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 24px;

    @include down(lg) {
      grid-template-columns: repeat(2, 1fr);
    }

    @include down(sm) {
      grid-template-columns: 1fr;
    }
  }
}

.tp-feature-card {
  background: var(--glass-background);
  border: 1px solid var(--glass-border);
  border-radius: 20px;
  padding: 32px 28px;
  backdrop-filter: blur(var(--glass-blur));
  -webkit-backdrop-filter: blur(var(--glass-blur));
  transition: transform 0.2s ease, box-shadow 0.2s ease;

  &:hover {
    transform: translateY(-4px);
    box-shadow: var(--glass-shadow-elevated);
  }

  &__icon {
    width: 52px;
    height: 52px;
    display: flex;
    align-items: center;
    justify-content: center;
    border-radius: 14px;
    background: var(--surface-brand-extra-light);
    color: var(--surface-brand-regular);
    margin-bottom: 20px;
  }

  h4 {
    margin-bottom: 8px;
  }

  p {
    color: var(--text-medium);
    line-height: 160%;
  }
}

// ─── HOW IT WORKS ──────────────────────────────────
.tp-how-it-works {
  width: 100%;
}

.tp-steps {
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 32px;
  position: relative;

  @include down(lg) {
    grid-template-columns: repeat(2, 1fr);
  }

  @include down(sm) {
    grid-template-columns: 1fr;
  }
}

.tp-step {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  position: relative;

  // connector line between steps (desktop only)
  &:not(:last-child)::after {
    content: '';
    position: absolute;
    top: 28px;
    left: calc(50% + 36px);
    width: calc(100% - 72px + 32px); // card width minus circle + gap
    height: 2px;
    background: linear-gradient(90deg, var(--surface-brand-regular), var(--surface-brand-light));
    pointer-events: none;

    @include down(lg) {
      display: none;
    }
  }

  &__number {
    width: 56px;
    height: 56px;
    border-radius: 50%;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 22px;
    font-weight: 900;
    color: var(--text-inverse);
    background: var(--surface-brand-regular);
    margin-bottom: 20px;
    position: relative;
    z-index: 1;
    box-shadow: 0 4px 20px color-mix(in srgb, var(--surface-brand-regular) 40%, transparent);
  }

  &__content {
    h4 {
      margin-bottom: 8px;
    }

    p {
      color: var(--text-medium);
      line-height: 160%;
      max-width: 260px;
    }
  }
}

// ─── STATS BANNER ──────────────────────────────────
.tp-stats-banner {
  width: 100%;
  display: flex;
  justify-content: center;
  gap: 16px;
  flex-wrap: wrap;
  background: var(--surface-brand-regular);
  border-radius: 24px;
  padding: 48px 40px;

  @include down(sm) {
    padding: 32px 24px;
  }
}

.tp-stat {
  flex: 1 1 180px;
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  gap: 4px;
  padding: 12px 8px;

  &__number {
    font-size: 40px;
    font-weight: 900;
    color: var(--surface-brand-extra-light);
    letter-spacing: -1.5px;
    line-height: 1.1;

    @include down(sm) {
      font-size: 32px;
    }
  }

  &__label {
    color: var(--surface-brand-light);
    white-space: nowrap;
  }
}

// ─── MODELS SHOWCASE ───────────────────────────────
.tp-models {
  width: 100%;

  &__grid {
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 32px;

    @include down(md) {
      grid-template-columns: 1fr;
    }
  }
}

.tp-model-card {
  background: var(--glass-background);
  border: 1px solid var(--glass-border);
  border-radius: 20px;
  padding: 36px 32px;
  backdrop-filter: blur(var(--glass-blur));
  -webkit-backdrop-filter: blur(var(--glass-blur));
  display: flex;
  flex-direction: column;
  gap: 24px;

  &__header {
    display: flex;
    flex-direction: column;
    gap: 8px;
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

  &__details {
    list-style: none;
    padding: 0;
    margin: 0;
    display: flex;
    flex-direction: column;
    gap: 12px;

    li {
      display: flex;
      align-items: center;
      gap: 10px;
      color: var(--text-medium);

      .tp-icon {
        color: var(--surface-brand-regular);
        flex-shrink: 0;
      }
    }
  }
}

// ─── CTA BANNER ────────────────────────────────────
.tp-cta-banner {
  width: 100%;
  border-radius: 24px;
  background:
    radial-gradient(ellipse at 20% 50%, color-mix(in srgb, var(--surface-brand-regular) 14%, transparent) 0%, transparent 70%),
    radial-gradient(ellipse at 80% 50%, color-mix(in srgb, var(--surface-brand-medium) 10%, transparent) 0%, transparent 70%),
    var(--glass-background);
  border: 1px solid var(--glass-border);
  backdrop-filter: blur(var(--glass-blur-strong));
  -webkit-backdrop-filter: blur(var(--glass-blur-strong));
  overflow: hidden;

  &__content {
    padding: 64px 48px;
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    gap: 16px;

    h2 {
      background: linear-gradient(135deg, var(--text-brand-regular), var(--text-brand-dark));
      background-clip: text;
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
    }

    p {
      max-width: 540px;
      color: var(--text-medium);
    }

    .row {
      margin-top: 8px;
    }

    @include down(sm) {
      padding: 40px 24px;
    }
  }
}

// ─── FAQ ───────────────────────────────────────────
.tp-faq {
  width: 100%;
  max-width: 780px;
  margin-inline: auto;
  padding-bottom: 120px;

  &__list {
    display: flex;
    flex-direction: column;
    gap: 8px;
  }

  &__item {
    border-radius: 12px;
    overflow: hidden;
    background: var(--glass-background);
    border: 1px solid var(--glass-border);
    backdrop-filter: blur(var(--glass-blur));
    -webkit-backdrop-filter: blur(var(--glass-blur));
    transition: box-shadow 0.2s ease;

    &--open {
      box-shadow: var(--glass-shadow);
    }
  }

  &__question {
    width: 100%;
    display: flex;
    justify-content: space-between;
    align-items: center;
    background: none;
    border: none;
    padding: 20px 24px;
    cursor: pointer;
    text-align: left;
    gap: 16px;
    color: var(--text);

    .tp-icon {
      flex-shrink: 0;
      color: var(--surface-brand-regular);
    }

    &:hover {
      background: var(--glass-background-light);
    }
  }

  &__answer {
    padding: 0 24px 20px;
    overflow: hidden;

    p {
      color: var(--text-medium);
      line-height: 170%;
    }
  }
}

// FAQ expand transition
.faq-expand-enter-active,
.faq-expand-leave-active {
  transition: all 0.25s ease;
}

.faq-expand-enter-from,
.faq-expand-leave-to {
  opacity: 0;
  transform: translateY(-8px);
}
</style>
