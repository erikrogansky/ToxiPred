<template>
  <div ref="container" class="w-full" />
</template>

<script setup lang="ts">
import { ref, watch, onMounted, nextTick } from 'vue'

const props = defineProps<{
  smiles: string
  atomScores: number[]
  title?: string
}>()

const container = ref<HTMLDivElement | null>(null)

async function render() {
  if (!container.value) return
  if (!props.smiles || !Array.isArray(props.atomScores) || props.atomScores.length === 0) return

  if (window.__xsmilesReady) await window.__xsmilesReady
  if (window.__rdkitReady) await window.__rdkitReady
  if (!window.RDKit || !window.xsmiles) return

  try {
    const expected = window.xsmiles.getNAtomsFromSmilesString(props.smiles)
    if (expected !== props.atomScores.length) {
      console.warn(`XSMILES expects ${expected} atom scores, got ${props.atomScores.length}`)
    }
  } catch {
    // ignore mismatch / failures
  }

  container.value.innerHTML = ''

  const setup = {
    molecule: {
      string: props.smiles,
      method: {
        name: 'atom_mask_dummy',
        scores: props.atomScores,
        attributes: {}
      },
      attributes: {}
    },
    width: 600,
    height: 300,
    bondLength: 75,
    gradientConfig: {
      palette: {
        name: 'MyPalette',
        colors: ['#7FFFD4', '#C7FFEA', '#FFFFFF', '#FCDDDD', '#F9B4B4']
      },
      colorDomain: [-0.09, -0.05, 0, 0.05, 0.09],
      highlight: false,
      blur: 0,
      opacity: { min: 0.1, max: 1.0 },
      radius: 16,
      thresholds: []
    },
    drawerType: 'RDKitDrawer',
    hideBarChart: false,
    hideAttributesTable: false,
    showScoresOnStructure: true
  }

  const api = window.xsmiles
  const fn = api.appendSingleView || api.appendSingle || api.renderSingle || api.appendView

  if (typeof fn !== 'function') {
    console.error('XSMILES API not found:', Object.keys(api))
    return
  }

  try {
    fn(container.value, setup)
  } catch (e) {
    console.error('XSMILES render failed. Setup was:', setup, e)
  }
}

onMounted(async () => {
  await nextTick()
  await render()
})

watch(
  () => [props.smiles, props.atomScores],
  async () => {
    await render()
  }
)
</script>

<style scoped>
div {
  min-height: 320px;
}
</style>
