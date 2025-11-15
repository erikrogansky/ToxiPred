<template>
  <div class="tp-feature-shap">
    <h2 class="tp-section-heading q-mb-md">
      Descriptor contributions
    </h2>

    <div v-if="!hasData" class="text-caption text-grey-6">
      No feature importance available for this prediction.
    </div>

    <div v-else ref="chartEl" class="tp-feature-shap__chart"></div>

    <div class="tp-feature-shap__legend text-caption text-grey-6 q-mt-sm">
      <span class="tp-dot tp-dot--positive"></span> Positive: increases predicted toxicity
      &nbsp;·&nbsp;
      <span class="tp-dot tp-dot--negative"></span> Negative: decreases predicted toxicity
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, computed, onMounted, onBeforeUnmount, watch } from 'vue'

type PlotlyLike = {
  newPlot: (
    el: HTMLElement,
    data: unknown[],
    layout?: Record<string, unknown>,
    config?: Record<string, unknown>,
  ) => unknown
  purge: (el: HTMLElement) => unknown
}

declare global {
  interface Window {
    Plotly?: PlotlyLike
    __plotlyReady?: Promise<void>
  }
}

const props = defineProps<{
  features: string[];
  scores: number[];
  prediction: number | null;
}>();

const chartEl = ref<HTMLDivElement | null>(null);

const hasData = computed(() => {
  const pred = props.prediction
  return (
    props.features?.length > 0 &&
    props.scores?.length > 0 &&
    pred !== null &&
    !Number.isNaN(pred)
  )
});

const sorted = computed(() => {
  const len = Math.min(props.features.length, props.scores.length);
  const pairs: { name: string; score: number }[] = [];
  for (let i = 0; i < len; i++) {
    const s = props.scores[i];
    if (s == null || Number.isNaN(s)) continue;
    pairs.push({ name: props.features[i] ?? '', score: s });
  }
  return pairs.sort((a, b) => Math.abs(b.score) - Math.abs(a.score));
});


async function draw() {
  if (!chartEl.value || !hasData.value) return;

  if (window.__plotlyReady) {
    await window.__plotlyReady;
  }

  const plotly = window.Plotly;
  if (!plotly) {
    console.error('Plotly not available');
    return;
  }

  const feats = sorted.value;

  if (!feats.length) return;

  const labels = feats.map(f => f.name);

  const measures = feats.map(() => 'relative' as const);

  const xValues = feats.map(f => f.score);

  const text = feats.map(f =>
    f.score >= 0 ? `+${f.score.toFixed(3)}` : f.score.toFixed(3),
  );

  const colors = feats.map(f =>
    f.score >= 0 ? '#ef4444' : '#3b82f6',
  );

  const data: unknown[] = [
    {
      type: 'waterfall',
      orientation: 'h',
      y: labels,
      x: xValues,
      measure: measures,
      text,
      textposition: 'outside',
      connector: {
        line: { color: 'rgba(0,0,0,0.25)' },
      },
      marker: { color: colors },
      decreasing: { marker: { color: '#3b82f6' } },
      increasing: { marker: { color: '#ef4444' } },
      hovertemplate: '%{y}<br>Δ: %{x:.4f}<extra></extra>',
    },
  ];

  const layout: Record<string, unknown> = {
    margin: { l: 140, r: 40, t: 10, b: 40 },
    xaxis: {
      title: 'Cumulative contribution',
      zeroline: true,
      zerolinecolor: 'rgba(0,0,0,0.2)',
      zerolinewidth: 1,
    },
    yaxis: {
      automargin: true,
    },
    showlegend: false,
    plot_bgcolor: 'rgba(0,0,0,0)',
    paper_bgcolor: 'rgba(0,0,0,0)',
  };

  const config: Record<string, unknown> = {
    displayModeBar: false,
    responsive: true,
  };

  await plotly.newPlot(chartEl.value, data, layout, config);
}


onMounted(() => {
  void draw();
});

watch(
  () => [props.features, props.scores, props.prediction],
  () => {
    void draw();
  },
  { deep: true },
);

onBeforeUnmount(() => {
  const plotly = window.Plotly;
  if (chartEl.value && plotly) {
    plotly.purge(chartEl.value);
  }
});
</script>

<style scoped lang="scss">
.tp-feature-shap {
  max-width: 960px;
  margin-inline: auto;
}

.tp-section-heading {
  font-size: 1.2rem;
  font-weight: 600;
}

.tp-feature-shap__chart {
  width: 100%;
  height: 420px;
}

.tp-feature-shap__legend {
  display: flex;
  justify-content: center;
  gap: 0.5rem;
  align-items: center;
  flex-wrap: wrap;
}

.tp-dot {
  display: inline-block;
  width: 8px;
  height: 8px;
  border-radius: 999px;
  margin-right: 4px;
}

.tp-dot--positive {
  background: #ef4444;
}

.tp-dot--negative {
  background: #3b82f6;
}
</style>
