<template>
  <q-dialog :model-value="modelValue ?? false" @update:model-value="emit('update:modelValue', $event)" persistent>
    <q-card class="tp-new-prediction-dialog">
      <q-card-section class="tp-dialog-header">
        <h3 class="tp-h3">New Prediction</h3>
        <tp-icon-button
          icon-name="close-square"
          weight="regular"
          :size="24"
          @click="close"
        />
      </q-card-section>

      <q-card-section class="tp-dialog-content">
        <!-- Tab Selection -->
        <tp-button-group 
          :labels="['Create New', 'Import from Link']" 
          :active-index="activeTab" 
          @click="activeTab = $event" 
          class="q-mb-lg"
        />

        <!-- Create New Tab -->
        <div v-if="activeTab === 0" class="tp-create-new-section">
          <p class="paragraph-regular q-mb-md">Enter a compound identifier to predict its toxicity.</p>
          
          <tp-input-with-select 
            select-label="Prediction" 
            label="Enter SMILES, CAS number, or trivial name" 
            hint="For batch submissions, separate the identifiers by a comma" 
            v-model:input-value="inputValue" 
            v-model:selected-value="selectedValue" 
            :options="modelOptions" 
            @enter="submitPrediction"
          />

          <div v-if="predictionError" class="tp-error-message q-mt-sm">
            {{ predictionError }}
          </div>

          <div class="tp-actions q-mt-lg">
            <tp-button
              label="Draw molecule"
              variant="outline"
              @click="goToDraw"
            />
            <div style="flex: 1" />
            <tp-button
              label="Cancel"
              variant="outline"
              @click="close"
            />
            <tp-button
              label="Predict"
              :disabled="!inputValue || !selectedValue"
              :loading="predictionLoading"
              @click="submitPrediction"
            />
          </div>
        </div>

        <!-- Import from Link Tab -->
        <div v-if="activeTab === 1" class="tp-import-section">
          <p class="paragraph-regular q-mb-md">Import a shared prediction into your workspace using a share link and password.</p>
          
          <div class="tp-input-group q-mb-md">
            <label class="tp-label">Share Link</label>
            <input
              v-model="shareLink"
              type="text"
              placeholder="Paste the full share link here"
              class="tp-text-input"
            />
            <span class="tp-hint">e.g., https://toxipred.com/shared/abc123xyz</span>
          </div>

          <div class="tp-input-group q-mb-md">
            <label class="tp-label">Password</label>
            <div class="tp-password-input-wrapper">
              <input
                v-model="importPassword"
                :type="showImportPassword ? 'text' : 'password'"
                placeholder="Enter password"
                class="tp-text-input"
                @keyup.enter="importFromLink"
              />
              <tp-icon-button
                :icon-name="showImportPassword ? 'eye-slash' : 'eye'"
                weight="regular"
                :size="18"
                @click="showImportPassword = !showImportPassword"
              />
            </div>
          </div>

          <div v-if="importError" class="tp-error-message q-mb-md">
            {{ importError }}
          </div>

          <div class="tp-actions">
            <tp-button
              label="Cancel"
              variant="outline"
              @click="close"
            />
            <tp-button
              label="Import to Workspace"
              :disabled="!shareLink || !importPassword"
              :loading="importLoading"
              @click="importFromLink"
            />
          </div>
        </div>
      </q-card-section>
    </q-card>
  </q-dialog>
</template>

<script setup lang="ts">
import { ref, computed } from 'vue';
import TpButton from 'components/TpButton.vue';
import TpButtonGroup from 'components/TpButtonGroup.vue';
import TpIconButton from 'components/TpIconButton.vue';
import TpInputWithSelect from 'components/TpInputWithSelect.vue';
import { useModelsStore, type TestType, type PredictionTarget } from 'src/stores/models-store';
import { useJobsStore } from 'src/stores/jobs-store';
import { api } from 'src/boot/axios';
import { useRouter } from 'vue-router';
const emit = defineEmits<{
  (e: 'update:modelValue', value: boolean): void;
}>();

defineProps<{
  modelValue?: boolean;
}>();

const modelsStore = useModelsStore();
const jobsStore = useJobsStore();
const router = useRouter();

const activeTab = ref(0);

// Create New tab state
const inputValue = ref('');
const selectedValue = ref('');
const predictionLoading = ref(false);
const predictionError = ref<string | null>(null);

// Import tab state
const shareLink = ref('');
const importPassword = ref('');
const showImportPassword = ref(false);
const importLoading = ref(false);
const importError = ref<string | null>(null);

// Extract token from share link
function extractTokenFromLink(link: string): string | null {
  // Handle various formats:
  // - Full URL: https://example.com/shared/abc123
  // - Path only: /shared/abc123
  // - Just token: abc123
  const trimmed = link.trim();
  
  // Try to match /shared/{token} pattern
  const match = trimmed.match(/\/shared\/([^/?#]+)/);
  if (match && match[1]) {
    return match[1];
  }
  
  // If no match and no slashes, assume it's just the token
  if (!trimmed.includes('/')) {
    return trimmed;
  }
  
  return null;
}

// Labels for display
const testTypeLabels: Record<TestType, string> = {
  'in_vitro': 'In Vitro',
  'in_vivo': 'In Vivo',
  'in_chemico': 'In Chemico',
};

const predictionTargetLabels: Record<PredictionTarget, string> = {
  'photo_irritation': 'Photo Irritation',
  'photo_toxicity': 'Phototoxicity',
};

const modelOptions = computed(() => {
  return modelsStore.getModels.map(modelName => {
    const detail = modelsStore.getModelDetail(modelName);
    const parts: string[] = [];
    
    if (detail?.prediction_target) {
      parts.push(predictionTargetLabels[detail.prediction_target] || detail.prediction_target);
    }
    if (detail?.test_type) {
      parts.push(testTypeLabels[detail.test_type] || detail.test_type);
    }
    
    const label = parts.length > 0 ? parts.join(' - ') : modelName;
    
    return {
      label,
      value: modelName
    };
  });
});

function goToDraw() {
  close();
  void router.push({ name: 'draw', query: { model: selectedValue.value || undefined } });
}

function close() {
  emit('update:modelValue', false);
  // Reset state
  activeTab.value = 0;
  inputValue.value = '';
  selectedValue.value = '';
  predictionError.value = null;
  shareLink.value = '';
  importPassword.value = '';
  importError.value = null;
}

async function submitPrediction() {
  if (!inputValue.value || !selectedValue.value) return;

  predictionLoading.value = true;
  predictionError.value = null;

  try {
    const response = await api.post(`/jobs/predict/${selectedValue.value}`, {
      query: inputValue.value
    });

    if (response.data.error) {
      predictionError.value = response.data.error;
      return;
    }

    if (!response.data.job_id) {
      predictionError.value = 'No job ID returned from server';
      return;
    }

    jobsStore.upsert({
      id: response.data.job_id,
      model: selectedValue.value,
      state: 'PENDING',
      percent: 0,
    });
    jobsStore.startTracking(response.data.job_id);

    close();
    // Stay on workspace page - the new job will appear in the list

  } catch (error) {
    console.error('Error making prediction:', error);
    predictionError.value = 'Error submitting prediction. Please try again.';
  } finally {
    predictionLoading.value = false;
  }
}

async function importFromLink() {
  if (!shareLink.value || !importPassword.value) return;

  const token = extractTokenFromLink(shareLink.value);
  if (!token) {
    importError.value = 'Invalid share link format. Please paste the full link.';
    return;
  }

  importLoading.value = true;
  importError.value = null;

  try {
    // First get the job data
    const dataRes = await api.post(`/jobs/shared/${token}`, {
      password: importPassword.value,
    });

    // Then import it
    const importRes = await api.post(`/jobs/import/${token}`, {
      password: importPassword.value,
    });

    const jobId = importRes.data.job_id;
    const jobData = dataRes.data.job_data;

    // Add to jobs store
    jobsStore.upsert({
      id: jobId,
      model: jobData?.model ?? '',
      name: jobData?.name ?? null,
      trivial_name: jobData?.trivial_name ?? null,
      formula: jobData?.formula ?? null,
      state: 'SUCCESS',
      prediction: Array.isArray(jobData?.prediction) 
        ? (jobData.prediction[0] as 0 | 1 | null) 
        : (jobData?.prediction as 0 | 1 | null) ?? null,
      canonical_smiles: jobData?.canonical_smiles ?? null,
      other_names: jobData?.other_names ?? null,
    });

    close();
    // Stay on workspace page - the imported job will appear in the list

  } catch (err: unknown) {
    const error = err as { response?: { status?: number; data?: { detail?: string } } };
    if (error.response?.status === 401) {
      importError.value = 'Incorrect password';
    } else if (error.response?.status === 404) {
      importError.value = 'Share link not found or has been deleted';
    } else if (error.response?.status === 410) {
      importError.value = error.response?.data?.detail || 'This prediction is no longer available';
    } else {
      importError.value = 'Failed to import prediction. Please check the token and password.';
    }
  } finally {
    importLoading.value = false;
  }
}
</script>

<style scoped lang="scss">
.tp-new-prediction-dialog {
  min-width: 480px;
  max-width: 600px;
  border-radius: 20px;
  background: var(--glass-background);
  backdrop-filter: var(--glass-blur);
  border: 1px solid var(--glass-border);
  box-shadow: var(--glass-shadow);
  
  @media (max-width: 600px) {
    min-width: unset;
    width: 95vw;
  }
}

.tp-dialog-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 20px 24px;
}

.tp-dialog-content {
  padding: 0 24px 24px;
}

.tp-input-group {
  display: flex;
  flex-direction: column;
  gap: 6px;
}

.tp-label {
  font-size: 14px;
  font-weight: 500;
  color: var(--text);
}

.tp-text-input {
  width: 100%;
  padding: 12px 16px;
  border: 1px solid var(--glass-border);
  border-radius: 10px;
  font-size: 14px;
  outline: none;
  background: var(--glass-background);
  backdrop-filter: blur(8px);
  transition: all 0.2s ease;

  &:focus {
    border-color: var(--stroke-brand-regular);
    box-shadow: 0 0 0 3px rgba(var(--color-brand-rgb), 0.1);
  }

  &::placeholder {
    color: var(--text-medium);
  }
}

.tp-hint {
  font-size: 12px;
  color: var(--text-medium);
}

.tp-password-input-wrapper {
  display: flex;
  align-items: center;
  gap: 8px;
  border: 1px solid var(--glass-border);
  border-radius: 10px;
  padding-right: 8px;
  background: var(--glass-background);
  backdrop-filter: blur(8px);
  transition: all 0.2s ease;
  
  &:focus-within {
    border-color: var(--stroke-brand-regular);
    box-shadow: 0 0 0 3px rgba(var(--color-brand-rgb), 0.1);
  }

  .tp-text-input {
    border: none;
    flex: 1;
    background: transparent;
    backdrop-filter: none;
  }
}

.tp-error-message {
  color: var(--text-negative);
  font-size: 14px;
}

.tp-actions {
  display: flex;
  justify-content: flex-end;
  gap: 12px;
}
</style>
