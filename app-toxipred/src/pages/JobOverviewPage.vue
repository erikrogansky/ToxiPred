<template>
  <tp-page class="column items-center text-center justify-center q-pa-xl">
    <div v-if="loading">
      Loading prediction details…
    </div>

    <div v-else-if="error">
      {{ error }}
    </div>

    <div v-else-if="result" class="column items-center">
      <h1 class="tp-main-heading q-mb-sm">
        {{ result.name || 'Unnamed compound' }}
      </h1>

      <div v-if="result.formula" class="q-mb-sm">
        <span v-html="formatChemFormula(result.formula)"></span>
      </div>

      <div class="q-mb-md text-caption">
        <span v-if="predictionValue !== null">
          Prediction:
          <strong>{{ predictionValue === 1 ? 'Toxic' : 'Non-toxic' }}</strong>
        </span>
        <span v-if="result.model">
          · Model: {{ result.model }}
        </span>
      </div>

      <tp-xsmiles-renderer
        v-if="smiles"
        :smiles="smiles"
        :atom-scores="atomScores"
        title="Atomic contributions"
      />

      <tp-feature-importance-renderer
        v-if="featureNames.length && featureScores.length && predictionValue !== null"
        class="q-mt-xl full-width"
        :features="featureNames"
        :scores="featureScores"
        :prediction="predictionValue"
      />

      <div v-if="!smiles" class="text-caption q-mt-md">
        No structure available for visualization.
      </div>
    </div>

    <div v-else>
      No response available.
    </div>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpXsmilesRenderer from 'components/TpXsmilesRenderer.vue';
import TpFeatureImportanceRenderer from 'components/TpFeatureImportanceRenderer.vue';
import { api } from 'src/boot/axios';
import { useRoute } from 'vue-router';
import { ref, onMounted, computed } from 'vue';

const route = useRoute();
const jobId = route.params.job_id?.toString() || '';

interface JobResultPayload {
  input_query?: string | null;
  input_type?: string | null;
  name?: string | null;
  formula?: string | null;
  canonical_smiles?: string | null;
  prediction?: number | number[] | null;
  confidence?: number[] | null;
  features_used?: string[] | null;
  model?: string | null;
  feature_scores?: number[] | null;
  atom_scores?: number[] | null;
}

const loading = ref(true);
const error = ref<string | null>(null);
const result = ref<JobResultPayload | null>(null);
const smiles = ref<string>('');
const atomScores = ref<number[]>([]);
const featureNames = ref<string[]>([]);
const featureScores = ref<number[]>([]);

onMounted(async () => {
  try {
    const response = await api.get(`/jobs/result/${jobId}`);
    if (response.data?.error) {
      error.value = response.data.error || 'Job is not finished or failed.';
      return;
    }

    result.value = response.data as JobResultPayload;
    smiles.value = result.value.canonical_smiles ?? '';
    atomScores.value = result.value.atom_scores ?? [];

    featureNames.value = result.value.features_used ?? [];
    featureScores.value = (result.value.feature_scores ?? []);

    console.log('Fetched job result:', result.value);
  } catch (err) {
    console.error('Error fetching job result:', err);
    error.value = 'Error fetching job result.';
  } finally {
    loading.value = false;
  }
});

const predictionValue = computed<number | null>(() => {
  if (!result.value || result.value.prediction == null) return null;
  const p = result.value.prediction;
  return Array.isArray(p) ? (p[0] ?? null) : p;
});

const subscriptMap: Record<string, string> = {
  '0': '₀',
  '1': '₁',
  '2': '₂',
  '3': '₃',
  '4': '₄',
  '5': '₅',
  '6': '₆',
  '7': '₇',
  '8': '₈',
  '9': '₉'
};

const formatChemFormula = (formula: string) => {
  if (!formula) return '';
  return formula.replace(/\d/g, digit => subscriptMap[digit] ?? digit);
};
</script>

<style scoped lang="scss">
</style>
