<template>
  <tp-page class="tp-home-page">
    <!-- <tp-button label="Click Me"></tp-button>
    <tp-button label="Click Me" variant="outline"></tp-button>
    <tp-button label="Click Me" variant="link"></tp-button> -->
    <h1 class="tp-main-heading">Predict dermatological toxicity from structure</h1>
    <h2 class="tp-h3 bold">OECD-aligned QSAR models with applicability domain & confidence.</h2>
    <tp-input-with-select v-model:input-value="inputValue" v-model:selected-value="selectedValue" :options="modelOptions" @enter="onEnter"/>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpInputWithSelect from 'components/TpInputWithSelect.vue';
import { ref } from 'vue';

import { useModelsStore } from 'src/stores/models-store';
import { useJobsStore } from 'src/stores/jobs-store';
import { storeToRefs } from 'pinia';

const modelsStore = useModelsStore()
const { getModels/*, isLoading */ } = storeToRefs(modelsStore)

const modelOptions = getModels.value.map(model => ({
  label: model,
  value: model
}))

const inputValue = ref("");
const selectedValue = ref("");

import { api } from 'src/boot/axios';
import { useRouter } from 'vue-router';
const router = useRouter();

const onEnter = async () => {
  try {
    const response = await api.post(`/jobs/predict/${selectedValue.value}`, {
      query: inputValue.value
    });

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
  }
};
</script>

<style scoped lang="scss">
.tp-home-page {
  padding-block: 128px;

  .tp-main-heading {
    margin-bottom: 16px;
  }
}
</style>
