<template>
  <tp-page class="tp-docs">
    <h1>Documentation</h1>
    <p class="paragraph tp-docs__subtitle">
      Learn how ToxiPred works — from submitting your first prediction to interpreting the results.
    </p>

    <nav class="tp-docs__tabs" role="tablist">
      <button
        v-for="tab in tabs"
        :key="tab.id"
        role="tab"
        :aria-selected="activeTab === tab.id"
        class="tp-docs__tab"
        :class="{ 'tp-docs__tab--active': activeTab === tab.id }"
        @click="activeTab = tab.id"
      >
        <tp-icon :icon-name="tab.icon" weight="regular" :size="18" />
        {{ tab.label }}
      </button>
    </nav>

    <Transition name="docs-fade" mode="out-in">
      <div :key="activeTab" class="tp-docs__content">

        <!-- Overview -->
        <template v-if="activeTab === 'overview'">
          <section>
            <h2>What is ToxiPred?</h2>
            <p class="paragraph">
              ToxiPred is a <strong>free, open-access research platform</strong> for predicting the
              dermatological toxicity of chemical compounds. It uses QSAR (Quantitative Structure–Activity
              Relationship) machine learning models developed and validated at the Slovak University of
              Technology in Bratislava (STU).
            </p>
            <p class="paragraph">
              The platform enables researchers, students, and regulatory professionals to quickly screen
              chemical compounds for potential skin-related toxicity — without any registration, software
              installation, or cost.
            </p>
            <div class="tp-docs__image-slot">
              <img :src="screenshotSources.overview" alt="ToxiPred application overview" />
            </div>
          </section>

          <section>
            <h2>Key Features</h2>
            <ul>
              <li><strong>Instant predictions</strong> — results typically arrive in under 5 seconds.</li>
              <li><strong>No registration required</strong> — use the platform immediately with zero setup.</li>
              <li><strong>Batch processing</strong> — submit multiple compounds at once and get results for all of them.</li>
              <li><strong>Explainable AI</strong> — every prediction includes SHAP-based feature importance scores and atom-level contribution maps.</li>
              <li><strong>Applicability domain</strong> — each prediction is flagged with whether the compound falls within the model's reliable prediction space.</li>
              <li><strong>Shareable results</strong> — generate time-limited, optionally password-protected links to share predictions with colleagues.</li>
              <li><strong>Privacy-first</strong> — no cookies, no analytics, no tracking. Your data is automatically deleted after 14 days.</li>
            </ul>
          </section>

          <section>
            <h2>Available Models</h2>
            <p class="paragraph">ToxiPred currently offers eight QSAR models across three dermatological toxicity endpoints:</p>
            <ul>
              <li><strong>Corrosion</strong> — in vitro 3D and in vivo models.</li>
              <li><strong>Skin irritation</strong> — in vitro 3D and in vivo rabbit models.</li>
              <li><strong>Phototoxicity</strong> — in chemico, in vitro 3T3, in vitro 3D, and in vivo models.</li>
            </ul>
            <p class="paragraph">
              Open the Models tab for dataset sizes, selected descriptors, fingerprints, validation metrics,
              thresholds, and hyperparameters for each model.
            </p>
          </section>
        </template>

        <!-- Getting Started -->
        <template v-else-if="activeTab === 'getting-started'">
          <section>
            <h2>Submitting a Prediction</h2>
            <p class="paragraph">
              You can submit a prediction from the <router-link to="/">Home page</router-link> or from your
              <router-link to="/workspace">Workspace</router-link>. There are three ways to identify a compound:
            </p>
            <ul>
              <li><strong>SMILES string</strong> — the most precise method. Enter a valid SMILES notation directly (e.g. <code>CC(=O)Oc1ccccc1C(=O)O</code> for Aspirin).</li>
              <li><strong>CAS number</strong> — enter a CAS registry number (e.g. <code>50-78-2</code>) and ToxiPred will resolve it automatically.</li>
              <li><strong>Compound name</strong> — enter a trivial or IUPAC name (e.g. "Aspirin") and it will be resolved to its SMILES representation.</li>
            </ul>
            <div class="tp-docs__image-slot">
              <img :src="screenshotSources.predictionInput" alt="Prediction input with model selector" />
            </div>
          </section>

          <section>
            <h2>Choosing a Model</h2>
            <p class="paragraph">
              Before submitting, select the toxicity endpoint you want to predict using the model selector
              dropdown. Each model targets a different toxicological endpoint and may produce different results
              for the same compound.
            </p>
          </section>

          <section>
            <h2>Drawing a Molecule</h2>
            <p class="paragraph">
              If you prefer to draw your compound, navigate to the <router-link to="/draw">Draw</router-link> page
              which provides an integrated Ketcher molecular editor. Draw or modify a structure visually, then
              submit it directly for prediction.
            </p>
            <div class="tp-docs__image-slot">
              <img :src="screenshotSources.ketcher" alt="Ketcher molecular editor" />
            </div>
          </section>

          <section>
            <h2>Batch Predictions</h2>
            <p class="paragraph">
              You can submit multiple compounds in a single job by entering one identifier per line (or
              separated by commas). All compounds will be processed together, and you can view results
              for each compound individually in the job overview.
            </p>
          </section>
        </template>

        <!-- Models -->
        <template v-else-if="activeTab === 'models'">
          <section>
            <h2>Model Documentation</h2>
            <p class="paragraph">
              These summaries are sourced from the final model deliverables. Runtime feature lists are treated
              as authoritative when a document and model artifact disagree.
            </p>
          </section>

          <section
            v-for="model in modelDocumentationList"
            :id="model.slug"
            :key="model.modelName"
            class="tp-docs__model-card"
          >
            <div class="tp-docs__model-header">
              <span class="tp-docs__model-badge">{{ model.badge }}</span>
              <div>
                <h3>{{ model.name }}</h3>
                <p class="paragraph-small">{{ model.algorithm }} for {{ model.endpoint.toLowerCase() }} prediction.</p>
              </div>
            </div>

            <dl class="tp-docs__model-facts">
              <div>
                <dt>Dataset</dt>
                <dd>{{ model.dataset }} ({{ model.observations }} compounds)</dd>
              </div>
              <div>
                <dt>Assay type</dt>
                <dd>{{ model.assayType }}</dd>
              </div>
              <div>
                <dt>Descriptors</dt>
                <dd>{{ model.descriptors }}</dd>
              </div>
              <div>
                <dt>Fingerprint</dt>
                <dd>{{ model.fingerprint }}</dd>
              </div>
              <div>
                <dt>Feature selection</dt>
                <dd>{{ model.featureSelection }}</dd>
              </div>
              <div>
                <dt>Threshold</dt>
                <dd>{{ model.threshold }}</dd>
              </div>
              <div>
                <dt>Tuning</dt>
                <dd>{{ model.tuning }}</dd>
              </div>
              <div>
                <dt>Class balance</dt>
                <dd>{{ model.classBalance.join('; ') }}</dd>
              </div>
            </dl>

            <h4>Validation Metrics</h4>
            <table class="tp-docs__metrics-table">
              <thead>
                <tr>
                  <th>Metric</th>
                  <th v-if="usesTrainTestMetrics(model)">Train</th>
                  <th v-if="usesTrainTestMetrics(model)">Test</th>
                  <th v-if="!usesTrainTestMetrics(model)">Stacking</th>
                  <th v-if="!usesTrainTestMetrics(model)">Voting</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="metric in model.metrics" :key="metric.label">
                  <td>{{ metric.label }}</td>
                  <td v-if="usesTrainTestMetrics(model)">{{ metric.train }}</td>
                  <td v-if="usesTrainTestMetrics(model)">{{ metric.test }}</td>
                  <td v-if="!usesTrainTestMetrics(model)">{{ metric.stacking }}</td>
                  <td v-if="!usesTrainTestMetrics(model)">{{ metric.voting }}</td>
                </tr>
              </tbody>
            </table>

            <h4>Hyperparameters</h4>
            <ul class="tp-docs__compact-list">
              <li v-for="item in model.hyperparameters" :key="item">{{ item }}</li>
            </ul>

            <div v-if="model.note" class="tp-docs__highlight">
              <p class="paragraph">{{ model.note }}</p>
            </div>
          </section>
        </template>

        <!-- Results -->
        <template v-else-if="activeTab === 'results'">
          <section>
            <h2>Understanding Your Results</h2>
            <p class="paragraph">
              After submitting a prediction, you'll be taken to the <strong>Job Overview</strong> page where
              you can see results for each compound. Every prediction includes:
            </p>
          </section>

          <section>
            <h3>Prediction Outcome</h3>
            <p class="paragraph">
              Each compound is classified as either <strong>toxic</strong> or <strong>non-toxic</strong> for the
              selected endpoint. The result is displayed prominently with a color-coded badge.
            </p>
          </section>

          <section>
            <h3>Confidence Score</h3>
            <p class="paragraph">
              The model's confidence in its prediction is shown as a percentage. Higher confidence means the
              model is more certain about the classification. Values near 50% indicate borderline cases where
              you should exercise additional caution.
            </p>
          </section>

          <section>
            <h3>Applicability Domain</h3>
            <p class="paragraph">
              ToxiPred checks whether your compound falls within the <strong>applicability domain</strong> of
              the model — the chemical space where the model's predictions are considered reliable. If a
              compound is outside the domain, the prediction is still provided but flagged with a warning.
            </p>
            <div class="tp-docs__highlight">
              <p class="paragraph">
                Predictions for compounds <strong>outside the applicability domain</strong> should be treated
                with extra caution, as the model has not seen similar structures during training.
              </p>
            </div>
          </section>

          <section>
            <h3>Feature Importance (SHAP)</h3>
            <p class="paragraph">
              Every prediction comes with a SHAP (SHapley Additive exPlanations) analysis showing which
              molecular descriptors contributed most to the prediction. This helps you understand
              <em>why</em> the model made its decision, not just <em>what</em> the decision was.
            </p>
            <div class="tp-docs__image-slot">
              <img :src="screenshotSources.featureImportance" alt="Descriptor contribution chart" />
            </div>
          </section>

          <section>
            <h3>Atom Contributions (XSMILES)</h3>
            <p class="paragraph">
              For a more granular view, ToxiPred provides atom-level contribution maps using the XSMILES
              renderer. Each atom in the molecule is color-coded to indicate whether it contributes positively
              or negatively to the predicted toxicity.
            </p>
            <div class="tp-docs__image-slot">
              <img :src="screenshotSources.atomContributions" alt="XSMILES atom contribution visualization" />
            </div>
          </section>

          <section>
            <h3>Molecular Descriptors</h3>
            <p class="paragraph">
              You can also inspect the computed molecular descriptors used by the model. These 32 descriptors
              capture various physicochemical and structural properties of the compound.
            </p>
          </section>
        </template>

        <!-- Workspace -->
        <template v-else-if="activeTab === 'workspace'">
          <section>
            <h2>Your Workspace</h2>
            <p class="paragraph">
              The <router-link to="/workspace">Workspace</router-link> is your personal dashboard for managing
              prediction jobs. All your submitted jobs are listed here with their status, creation time, and
              compound count.
            </p>
            <div class="tp-docs__image-slot">
              <img :src="screenshotSources.workspace" alt="Workspace overview" />
            </div>
          </section>

          <section>
            <h2>Job History</h2>
            <p class="paragraph">
              Your job history is stored locally in your browser using localStorage. This means:
            </p>
            <ul>
              <li>Your history is <strong>private</strong> — it never leaves your device.</li>
              <li>Clearing your browser data will remove your local history.</li>
              <li>Server-side results are automatically deleted after <strong>14 days</strong>.</li>
              <li>Jobs older than 14 days will still appear in your local history but their detailed results will no longer be available from the server.</li>
            </ul>
          </section>

          <section>
            <h2>Sharing Results</h2>
            <p class="paragraph">
              You can share any prediction job with others by generating a share link. When sharing:
            </p>
            <ul>
              <li>A unique, cryptographically secure link is created.</li>
              <li>You can optionally protect the link with a password.</li>
              <li>Share links expire automatically after <strong>30 days</strong>.</li>
              <li>Recipients do not need an account — they can view results directly via the link.</li>
            </ul>
          </section>

          <section>
            <h2>Filtering &amp; Sorting</h2>
            <p class="paragraph">
              Use the filter controls in the Workspace to narrow down results by prediction outcome, model
              type, applicability domain status, or search by compound name. You can also sort jobs by
              date or compound count.
            </p>
          </section>
        </template>

        <!-- FAQ -->
        <template v-else-if="activeTab === 'faq'">
          <section>
            <h2>Frequently Asked Questions</h2>
          </section>

          <section>
            <h3>Is ToxiPred free to use?</h3>
            <p class="paragraph">
              Yes, completely. ToxiPred is a free, open-access academic research tool. No registration,
              subscription, or payment is required.
            </p>
          </section>

          <section>
            <h3>Can I use ToxiPred results in regulatory submissions?</h3>
            <p class="paragraph">
              ToxiPred is intended for <strong>research and screening purposes only</strong>. While the
              models follow OECD validation principles, the results should not be used as the sole basis
              for regulatory decisions. Always consult with qualified toxicologists and refer to official
              regulatory guidelines.
            </p>
          </section>

          <section>
            <h3>What happens to my data?</h3>
            <p class="paragraph">
              Chemical data you submit is processed to generate predictions and cached on the server for
              up to 14 days. After that, it is automatically and permanently deleted. We do not collect
              any personal information, use cookies, or track users. See our
              <router-link to="/privacy-policy">Privacy Policy</router-link> for full details.
            </p>
          </section>

          <section>
            <h3>What is the applicability domain?</h3>
            <p class="paragraph">
              The applicability domain defines the chemical space where a QSAR model can make reliable
              predictions. It is determined by the structural diversity of the training data. Compounds
              that fall outside this domain may receive less accurate predictions, which is why ToxiPred
              flags them with a warning.
            </p>
          </section>

          <section>
            <h3>What SMILES format is accepted?</h3>
            <p class="paragraph">
              ToxiPred accepts standard SMILES notation. The input is parsed and canonicalized using RDKit,
              so most valid SMILES strings will work. Stereochemistry notation (e.g. <code>@@</code>) and
              aromatic representations are supported.
            </p>
          </section>

          <section>
            <h3>How accurate are the predictions?</h3>
            <p class="paragraph">
              Model performance varies by endpoint. Detailed validation metrics (accuracy, sensitivity,
              specificity, AUC) are available on the model cards shown on the Home page. In general,
              predictions for compounds within the applicability domain are more reliable than those outside it.
            </p>
          </section>

          <section>
            <h3>Can I submit proprietary compounds?</h3>
            <p class="paragraph">
              You can, but please review our <router-link to="/security">Security</router-link> page first.
              While we encrypt all data in transit and delete it automatically, ToxiPred is an academic
              project without a dedicated security team. For highly sensitive structures, consider using
              SMILES notation (which doesn't reveal the compound's identity without specialized knowledge).
            </p>
          </section>

          <section>
            <h3>Who developed ToxiPred?</h3>
            <p class="paragraph">
              ToxiPred was developed at the <strong>Slovak University of Technology in Bratislava (STU)</strong>,
              Faculty of Chemical and Food Technology. For questions or collaboration inquiries, contact
              <tp-link href="mailto:roganskyerik@gmail.com" label="roganskyerik@gmail.com" underline="hover" /> or
              <tp-link href="mailto:marta.prnova@stuba.sk" label="marta.prnova@stuba.sk" underline="hover" />.
            </p>
          </section>
        </template>
      </div>
    </Transition>
  </tp-page>
</template>

<script setup lang="ts">
import { nextTick, onMounted, ref, watch } from 'vue';
import { useRoute } from 'vue-router';
import TpPage from 'components/TpPage.vue';
import TpIcon from 'components/TpIcon.vue';
import TpLink from 'components/TpLink.vue';
import {
  modelDocumentationList,
  type ModelDocumentation,
} from 'src/data/modelDocumentation';

import type { IconName } from 'src/css/icons/icon-names';

interface Tab {
  id: string;
  label: string;
  icon: IconName;
}

const tabs: Tab[] = [
  { id: 'overview', label: 'Overview', icon: 'book' },
  { id: 'getting-started', label: 'Getting Started', icon: 'play' },
  { id: 'models', label: 'Models', icon: 'document' },
  { id: 'results', label: 'Results', icon: 'chart' },
  { id: 'workspace', label: 'Workspace', icon: 'settings' },
  { id: 'faq', label: 'FAQ', icon: 'star' },
];

const route = useRoute();
const activeTab = ref('overview');

const screenshotSources = {
  overview: '/docs/screenshots/application-overview.png',
  predictionInput: '/docs/screenshots/prediction-input.png',
  ketcher: '/docs/screenshots/ketcher-editor.png',
  featureImportance: '/docs/screenshots/descriptor-contributions.png',
  atomContributions: '/docs/screenshots/xsmiles-atom-contributions.png',
  workspace: '/docs/screenshots/workspace-overview.png',
};

function usesTrainTestMetrics(model: ModelDocumentation) {
  return model.metrics.some((metric) => metric.train !== undefined || metric.test !== undefined);
}

async function syncTabWithHash() {
  if (!route.hash.startsWith('#model-')) return;

  activeTab.value = 'models';
  await nextTick();

  const target = document.getElementById(route.hash.slice(1));
  target?.scrollIntoView({ behavior: 'smooth', block: 'start' });
}

onMounted(() => {
  void syncTabWithHash();
});

watch(
  () => route.hash,
  () => {
    void syncTabWithHash();
  },
);
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins' as *;
@use 'src/css/helpers/glass' as glass;

.tp-docs {
  box-sizing: content-box;
  max-width: 920px;
  gap: 0;

  &__subtitle {
    color: var(--text-medium);
    margin-bottom: 32px;
    line-height: 1.7;
  }

  // ── Tabs ──
  &__tabs {
    display: flex;
    gap: 8px;
    margin-bottom: 40px;
    flex-wrap: wrap;
  }

  &__tab {
    @include glass.glass-chip(36px);
    display: inline-flex;
    align-items: center;
    gap: 6px;
    font-size: 14px !important;
    padding: 8px 16px !important;

    &--active {
      @include glass.glass-chip-active();
    }
  }

  // ── Content area ──
  &__content {
    section {
      margin-bottom: 32px;
    }

    h2 {
      margin-bottom: 12px;
      padding-top: 8px;
    }

    h3 {
      margin-top: 16px;
      margin-bottom: 8px;
    }

    p {
      margin-bottom: 12px;
      color: var(--text);
      line-height: 1.7;
    }

    ul {
      padding-left: 24px;
      margin-bottom: 12px;

      li {
        color: var(--text);
        line-height: 1.8;
        font-size: 16px;
      }
    }

    a {
      color: var(--text-brand-regular);
      text-decoration: none;

      &:hover {
        text-decoration: underline;
      }
    }

    code {
      background: var(--glass-background);
      border: 1px solid var(--glass-border);
      border-radius: 6px;
      padding: 2px 6px;
      font-size: 14px;
      font-family: 'SF Mono', 'Fira Code', monospace;
    }
  }

  // ── Image placeholders ──
  &__image-slot {
    @include glass.glass($radius: 12px);
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 200px;
    margin: 16px 0 24px;
    overflow: hidden;

    img {
      width: 100%;
      height: auto;
      display: block;
      border-radius: 12px;
    }
  }

  &__image-placeholder {
    color: var(--text-medium);
    font-size: 14px;
    font-style: italic;
  }

  &__model-card {
    @include glass.glass($radius: 12px);
    padding: 24px;
    scroll-margin-top: 96px;

    h4 {
      margin: 20px 0 8px;
    }
  }

  &__model-header {
    display: flex;
    align-items: flex-start;
    gap: 16px;
    margin-bottom: 20px;

    h3 {
      margin: 0 0 4px;
    }

    p {
      margin: 0;
      color: var(--text-medium);
    }
  }

  &__model-badge {
    display: inline-block;
    flex-shrink: 0;
    font-size: 11px;
    font-weight: 900;
    text-transform: uppercase;
    letter-spacing: 0.6px;
    padding: 4px 12px;
    border-radius: 360px;
    background: var(--surface-brand-extra-light);
    color: var(--text-brand-regular);
  }

  &__model-facts {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 14px 20px;
    margin: 0;

    div {
      min-width: 0;
    }

    dt {
      font-size: 12px;
      font-weight: 700;
      color: var(--text-medium);
      text-transform: uppercase;
      letter-spacing: 0.5px;
      margin-bottom: 4px;
    }

    dd {
      margin: 0;
      color: var(--text);
      line-height: 1.5;
    }
  }

  &__metrics-table {
    width: 100%;
    border-collapse: collapse;
    overflow: hidden;
    border-radius: 8px;

    th,
    td {
      text-align: left;
      padding: 10px 12px;
      border-bottom: 1px solid var(--glass-border);
    }

    th {
      background: rgba(128, 128, 128, 0.15);
      color: var(--text-brand-regular);
      font-size: 13px;
    }

    td {
      color: var(--text);
      font-size: 14px;
    }
  }

  &__compact-list {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 6px 18px;
    padding-left: 18px !important;

    li {
      font-size: 14px !important;
      line-height: 1.5 !important;
    }
  }

  // ── Info highlight box ──
  &__highlight {
    @include glass.glass-info-box();
    border-left: 4px solid var(--accent-blue-regular);
    padding: 20px 24px;
    margin: 8px 0 16px;

    p {
      margin-bottom: 0;
    }
  }
}

// ── Tab content transition ──
.docs-fade-enter-active,
.docs-fade-leave-active {
  transition: opacity 0.2s ease;
}

.docs-fade-enter-from,
.docs-fade-leave-to {
  opacity: 0;
}

// ── Responsive ──
@include down(sm) {
  .tp-docs__tabs {
    gap: 6px;
  }

  .tp-docs__tab {
    font-size: 13px !important;
    padding: 6px 12px !important;
  }

  .tp-docs__image-slot {
    min-height: 140px;
  }

  .tp-docs__model-header {
    flex-direction: column;
  }

  .tp-docs__model-facts,
  .tp-docs__compact-list {
    grid-template-columns: 1fr;
  }
}
</style>
