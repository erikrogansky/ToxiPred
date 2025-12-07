<template>
  <tp-page class="tp-shared-job column items-center text-center justify-center q-pa-md q-pa-lg-xl">
    <!-- Password Entry State -->
    <div v-if="state === 'password'" class="tp-password-panel column items-center">
      <tp-icon icon-name="lock" weight="regular" :size="48" class="q-mb-md" />
      <h1 class="tp-main-heading q-mb-sm">Protected Prediction</h1>
      <p class="paragraph-regular q-mb-lg">Enter the password to view this shared prediction.</p>
      
      <div class="tp-password-input-wrapper q-mb-md">
        <input
          v-model="password"
          :type="showPassword ? 'text' : 'password'"
          placeholder="Enter password"
          class="tp-password-input"
          @keyup.enter="accessShared"
        />
        <tp-icon-button
          :icon-name="showPassword ? 'eye-slash' : 'eye'"
          weight="regular"
          :size="20"
          @click="showPassword = !showPassword"
        />
      </div>
      
      <div v-if="passwordError" class="tp-error-message q-mb-md">
        {{ passwordError }}
      </div>
      
      <div class="row q-gutter-sm">
        <tp-button
          label="View Prediction"
          :loading="loading"
          @click="accessShared"
        />
        <tp-button
          label="Import to Workspace"
          variant="outline"
          :loading="importLoading"
          @click="importShared"
        />
      </div>
    </div>

    <!-- Loading State -->
    <div v-else-if="state === 'loading'" class="column items-center">
      <q-spinner-dots size="48px" color="primary" />
      <p class="q-mt-md">Loading prediction...</p>
    </div>

    <!-- Error/Deleted State -->
    <div v-else-if="state === 'error'" class="tp-error-panel column items-center">
      <tp-icon icon-name="warning-2" weight="regular" :size="48" class="q-mb-md text-negative" />
      <h1 class="tp-main-heading q-mb-sm">{{ errorTitle }}</h1>
      <p class="paragraph-regular q-mb-lg">{{ errorMessage }}</p>
      <tp-button
        label="Go to Home"
        @click="$router.push('/')"
      />
    </div>

    <!-- Success: Show Prediction Result -->
    <div v-else-if="state === 'success' && result" class="column items-center full-width">
      <div class="tp-header full-width">
        <div class="tp-header__left">
          <tp-icon-button
            icon-name="arrow-left-02"
            weight="regular"
            :size="28"
            @click="$router.push('/')"
          />
          <h1 class="tp-main-heading">{{ displayName }}</h1>
          <span class="tp-shared-badge">Shared</span>
        </div>
        <div class="tp-header__right">
          <tp-button
            label="Import to Workspace"
            variant="outline"
            :loading="importLoading"
            @click="importShared"
          />
        </div>
      </div>

      <div class="tp-content-row">
        <tp-compound-info
          class="tp-panel"
          :trivial-name="result.trivial_name"
          :name="result.name"
          :formula="result.formula"
          :other-names="result.other_names"
          :smiles="result.canonical_smiles"
          :input-type="result.input_type"
        />

        <div class="tp-molecule-preview tp-panel">
          <tp-button-group :labels="buttonLabels" :active-index="selectedIndex" @click="selectedIndex = $event" />
          <tp-molecule-renderer 
            v-if="smiles" 
            :smiles="smiles" 
            :mode="selectedIndex === 0 ? '2d' : '3d'"
          />
        </div>
      </div>

      <div class="tp-content-row">
        <tp-descriptors-summary
          class="tp-panel"
          :features="featureNames"
          :values="featureValues"
          :scores="featureScores"
          :prediction="predictionValue"
          :confidence="confidence"
          :positive-label="positiveLabel"
          :negative-label="negativeLabel"
        />

        <div class="tp-visualization-panel tp-panel">
          <tp-button-group 
            :labels="visualizationLabels" 
            :active-index="visualizationIndex" 
            @click="visualizationIndex = $event" 
          />
          
          <tp-xsmiles-renderer
            v-if="visualizationIndex === 0 && smiles"
            :smiles="smiles"
            :atom-scores="atomScores"
            title="Atomic contributions"
          />

          <tp-feature-importance-renderer
            v-if="visualizationIndex === 1 && featureNames.length && featureScores.length && predictionValue !== null"
            class="full-width"
            :features="featureNames"
            :scores="featureScores"
            :prediction="predictionValue"
          />
        </div>
      </div>
    </div>

    <!-- Import Success Notification -->
    <q-dialog v-model="showImportSuccess">
      <q-card class="tp-import-success-card">
        <q-card-section class="column items-center text-center">
          <tp-icon 
            :icon-name="alreadyInWorkspace ? 'info-circle' : 'tick-circle'" 
            weight="regular" 
            :size="48" 
            :class="alreadyInWorkspace ? 'text-warning' : 'text-positive'" 
            class="q-mb-md" 
          />
          <h3 class="tp-h3 q-mb-sm">
            {{ alreadyInWorkspace ? 'Already in Workspace' : 'Imported Successfully' }}
          </h3>
          <p class="paragraph-regular">
            {{ alreadyInWorkspace 
              ? 'This prediction is already in your workspace.' 
              : 'The prediction has been added to your workspace.' 
            }}
          </p>
        </q-card-section>
        <q-card-actions align="center">
          <tp-button label="Go to Workspace" @click="goToWorkspace" />
          <tp-button label="Stay Here" variant="outline" @click="showImportSuccess = false" />
        </q-card-actions>
      </q-card>
    </q-dialog>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpXsmilesRenderer from 'components/TpXsmilesRenderer.vue';
import TpFeatureImportanceRenderer from 'components/TpFeatureImportanceRenderer.vue';
import TpDescriptorsSummary from 'components/TpDescriptorsSummary.vue';
import TpButtonGroup from 'components/TpButtonGroup.vue';
import TpMoleculeRenderer from 'components/TpMoleculeRenderer.vue';
import TpButton from 'components/TpButton.vue';
import TpIconButton from 'components/TpIconButton.vue';
import TpIcon from 'components/TpIcon.vue';
import TpCompoundInfo from 'components/TpCompoundInfo.vue';
import { api } from 'src/boot/axios';
import { useRoute, useRouter } from 'vue-router';
import { ref, onMounted, computed } from 'vue';
import { useModelsStore } from 'src/stores/models-store';
import { useJobsStore } from 'src/stores/jobs-store';

const route = useRoute();
const router = useRouter();
const token = route.params.token?.toString() || '';
const modelsStore = useModelsStore();
const jobsStore = useJobsStore();

type ViewState = 'password' | 'loading' | 'success' | 'error';

interface JobResultPayload {
  input_query?: string | null;
  input_type?: string | null;
  name?: string | null;
  trivial_name?: string | null;
  formula?: string | null;
  canonical_smiles?: string | null;
  other_names?: string[] | null;
  prediction?: number | number[] | null;
  confidence?: number | null;
  features_used?: string[] | null;
  feature_values?: (number | null)[] | null;
  model?: string | null;
  feature_scores?: number[] | null;
  atom_scores?: number[] | null;
}

const state = ref<ViewState>('loading');
const password = ref('');
const showPassword = ref(false);
const passwordError = ref<string | null>(null);
const loading = ref(false);
const importLoading = ref(false);
const errorTitle = ref('');
const errorMessage = ref('');
const showImportSuccess = ref(false);

const result = ref<JobResultPayload | null>(null);
const smiles = ref<string>('');
const atomScores = ref<number[]>([]);
const featureNames = ref<string[]>([]);
const featureValues = ref<(number | null)[]>([]);
const featureScores = ref<number[]>([]);
const confidence = ref<number | null>(null);

onMounted(async () => {
  // Check if share link is valid before asking for password
  try {
    const checkRes = await api.get(`/jobs/shared/${token}/check`);
    if (!checkRes.data.valid) {
      handleInvalidShare(checkRes.data.reason);
      return;
    }
    state.value = 'password';
  } catch {
    handleError('Connection Error', 'Unable to verify share link. Please try again later.');
  }
});

function handleInvalidShare(reason: string) {
  switch (reason) {
    case 'not_found':
      handleError('Link Not Found', 'This share link does not exist or has been removed.');
      break;
    case 'expired':
      handleError('Link Expired', 'This share link has expired.');
      break;
    case 'deleted':
      handleError('Prediction Deleted', 'The shared prediction has been deleted by all users who had access to it.');
      break;
    case 'job_expired':
      handleError('Prediction Expired', 'The prediction data has expired and is no longer available.');
      break;
    default:
      handleError('Invalid Link', 'This share link is no longer valid.');
  }
}

function handleError(title: string, message: string) {
  errorTitle.value = title;
  errorMessage.value = message;
  state.value = 'error';
}

async function accessShared() {
  if (!password.value.trim()) {
    passwordError.value = 'Please enter the password';
    return;
  }

  loading.value = true;
  passwordError.value = null;

  try {
    const response = await api.post(`/jobs/shared/${token}`, {
      password: password.value,
    });

    result.value = response.data.job_data as JobResultPayload;
    smiles.value = result.value.canonical_smiles ?? '';
    atomScores.value = result.value.atom_scores ?? [];
    featureNames.value = result.value.features_used ?? [];
    featureValues.value = result.value.feature_values ?? [];
    featureScores.value = result.value.feature_scores ?? [];
    confidence.value = result.value.confidence ?? null;

    state.value = 'success';
  } catch (err: unknown) {
    const error = err as { response?: { status?: number; data?: { detail?: string } } };
    if (error.response?.status === 401) {
      passwordError.value = 'Incorrect password';
    } else if (error.response?.status === 410) {
      handleError('No Longer Available', error.response?.data?.detail || 'This prediction is no longer available.');
    } else {
      passwordError.value = 'Failed to access shared prediction';
    }
  } finally {
    loading.value = false;
  }
}

const alreadyInWorkspace = ref(false);

async function importShared() {
  if (!password.value.trim() && state.value === 'password') {
    passwordError.value = 'Please enter the password';
    return;
  }

  importLoading.value = true;
  passwordError.value = null;

  try {
    // Send existing job IDs so backend can check for duplicates
    const existingJobIds = Array.from(jobsStore.jobs.keys());
    
    const response = await api.post(`/jobs/import/${token}`, {
      password: password.value,
      existing_job_ids: existingJobIds,
    });

    const jobId = response.data.job_id;
    const alreadyExists = response.data.already_exists;
    
    if (alreadyExists) {
      // Job already in workspace - show different message
      alreadyInWorkspace.value = true;
      showImportSuccess.value = true;
      return;
    }
    
    // Also fetch the job data if we haven't already
    if (!result.value) {
      const dataRes = await api.post(`/jobs/shared/${token}`, {
        password: password.value,
      });
      result.value = dataRes.data.job_data as JobResultPayload;
    }

    // Add to jobs store
    jobsStore.upsert({
      id: jobId,
      model: result.value?.model ?? '',
      name: result.value?.name ?? null,
      trivial_name: result.value?.trivial_name ?? null,
      formula: result.value?.formula ?? null,
      state: 'SUCCESS',
      prediction: Array.isArray(result.value?.prediction) 
        ? (result.value.prediction[0] as 0 | 1 | null) 
        : (result.value?.prediction as 0 | 1 | null) ?? null,
      canonical_smiles: result.value?.canonical_smiles ?? null,
      other_names: result.value?.other_names ?? null,
    });

    alreadyInWorkspace.value = false;
    showImportSuccess.value = true;
  } catch (err: unknown) {
    const error = err as { response?: { status?: number; data?: { detail?: string } } };
    if (error.response?.status === 401) {
      passwordError.value = 'Incorrect password';
    } else if (error.response?.status === 410) {
      handleError('No Longer Available', error.response?.data?.detail || 'This prediction is no longer available.');
    } else {
      passwordError.value = 'Failed to import prediction';
    }
  } finally {
    importLoading.value = false;
  }
}

function goToWorkspace() {
  void router.push({ name: 'workspace' });
}

const predictionValue = computed<number | null>(() => {
  if (!result.value || result.value.prediction == null) return null;
  const p = result.value.prediction;
  return Array.isArray(p) ? (p[0] ?? null) : p;
});

const displayName = computed<string>(() => {
  if (!result.value) return 'Shared Prediction';
  return result.value.trivial_name || result.value.canonical_smiles || result.value.name || 'Unknown Compound';
});

const positiveLabel = computed<string>(() => {
  if (!result.value) return 'Positive';
  const modelName = result.value.model;
  if (modelName) {
    const detail = modelsStore.getModelDetail(modelName);
    return detail?.positive_label || 'Toxic';
  }
  return 'Toxic';
});

const negativeLabel = computed<string>(() => {
  if (!result.value) return 'Negative';
  const modelName = result.value.model;
  if (modelName) {
    const detail = modelsStore.getModelDetail(modelName);
    return detail?.negative_label || 'Non-toxic';
  }
  return 'Non-toxic';
});

const buttonLabels = ['2D Structure', '3D Structure'];
const selectedIndex = ref(0);

const visualizationLabels = ['XSMILES', 'Descriptor Importance'];
const visualizationIndex = ref(0);
</script>

<style scoped lang="scss">
.tp-shared-job {
  --content-gap: 32px;
  
  @media (min-width: 1024px) {
    --content-gap: 48px;
  }
}

.tp-password-panel {
  max-width: 400px;
  padding: 48px 32px;
  background: var(--surface-white);
  border-radius: 16px;
  border: 1px solid var(--stroke-extra-light);
}

.tp-password-input-wrapper {
  display: flex;
  align-items: center;
  gap: 8px;
  border: 1px solid var(--stroke-light);
  border-radius: 8px;
  padding: 8px 12px;
  width: 100%;
  max-width: 300px;
  
  &:focus-within {
    border-color: var(--stroke-brand-regular);
  }
}

.tp-password-input {
  flex: 1;
  border: none;
  outline: none;
  background: transparent;
  font-size: 14px;
}

.tp-error-message {
  color: var(--text-negative);
  font-size: 14px;
}

.tp-error-panel {
  max-width: 400px;
  padding: 48px 32px;
}

.tp-shared-badge {
  background: var(--surface-brand-extra-light);
  color: var(--text-brand-regular);
  padding: 4px 12px;
  border-radius: 16px;
  font-size: 12px;
  font-weight: 500;
}

.tp-header {
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  gap: 16px;
  margin-bottom: 28px;
  width: 100%;

  @media (min-width: 768px) {
    flex-direction: row;
    align-items: center;
    justify-content: space-between;
  }

  &__left {
    display: flex;
    align-items: center;
    gap: 16px;
  }

  &__right {
    display: flex;
    gap: 12px;
    flex-wrap: wrap;
  }
}

.tp-content-row {
  display: flex;
  flex-direction: column;
  align-items: stretch;
  width: 100%;
  gap: var(--content-gap);
  margin-bottom: var(--content-gap);

  @media (min-width: 900px) {
    flex-direction: row;
  }
}

.tp-panel {
  flex: 1;
  min-width: 0;
}

.tp-molecule-preview {
  display: flex;
  flex-direction: column;
  gap: 16px;
  border-radius: 12px;
  height: 100%;
}

.tp-visualization-panel {
  display: flex;
  flex-direction: column;
  gap: 16px;
  border-radius: 12px;
  height: 100%;
}

.tp-import-success-card {
  padding: 24px;
  min-width: 320px;
}
</style>
