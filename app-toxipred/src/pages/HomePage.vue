<template>
  <tp-page class="tp-home-page column items-center">
    <h1 class="tp-main-heading">Predict dermatological toxicity from structure</h1>
    <h2 class="tp-h3 bold">OECD-aligned QSAR models with applicability domain & confidence.</h2>

    <div class="row items-center justify-center full-width q-mt-md q-gutter-md">
      <tp-button label="See demo predictions" />
      <tp-button label="Documentation" variant="outline" />
    </div>

    <div class="row items-top justify-center full-width tp-input-section">
      <tp-input-with-select select-label="Prediction" label="Enter SMILES, CAS number, or trivial name" hint="For batch submissions, separate the identifiers by a comma" v-model:input-value="inputValue" v-model:selected-value="selectedValue" :options="modelOptions" @enter="submitPrediction"/>
      <tp-button :disabled="!inputValue || !selectedValue" style="margin: 0 0 20px 16px;" @click="submitPrediction" label="Predict" variant="outline" size="regular"></tp-button>
    </div>
</tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpInputWithSelect from 'components/TpInputWithSelect.vue';
import TpButton from 'components/TpButton.vue';
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

    // Check if there's an error in the response
    if (response.data.error) {
      console.error('Prediction error:', response.data.error);
      alert(`Error: ${response.data.error}`);
      return;
    }

    // Check if job_id exists
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
</script>

<style scoped lang="scss">
.tp-home-page {
  padding-block: 128px;

  .tp-main-heading {
    margin-bottom: 16px;
  }

  .tp-input-section {
    margin-top: 72px;
  }
}
</style>
