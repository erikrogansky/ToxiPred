<template>
  <div class="tp-molecule-renderer">
    <div v-if="mode === '2d'" ref="canvas2d" class="tp-molecule-renderer__canvas"></div>
    <div v-else-if="mode === '3d'" ref="canvas3d" class="tp-molecule-renderer__canvas tp-molecule-renderer__canvas--3d"></div>
    <div v-if="error" class="tp-molecule-renderer__error text-caption text-negative">
      {{ error }}
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, watch, onMounted, nextTick } from 'vue';

interface Props {
  smiles: string;
  mode: '2d' | '3d';
}

const props = defineProps<Props>();

const canvas2d = ref<HTMLDivElement | null>(null);
const canvas3d = ref<HTMLDivElement | null>(null);
const error = ref<string | null>(null);

// eslint-disable-next-line @typescript-eslint/no-explicit-any
let viewer3d: any = null;

const render2D = async () => {
  if (!canvas2d.value || !props.smiles) return;
  
  error.value = null;
  
  try {
    if (window.__rdkitReady) await window.__rdkitReady;
    if (!window.RDKit) {
      error.value = 'RDKit not available';
      return;
    }

    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const RDKit = window.RDKit as any;
    const mol = RDKit.get_mol(props.smiles);
    
    if (!mol || !mol.is_valid()) {
      error.value = 'Invalid SMILES structure';
      return;
    }

    // Get container width and calculate proportional height (4:3 aspect ratio)
    const containerWidth = canvas2d.value.clientWidth;
    const height = Math.round(containerWidth * 0.75);

    const svg = mol.get_svg(containerWidth, height);
    canvas2d.value.innerHTML = svg;
    
    mol.delete();
  } catch (err) {
    console.error('Error rendering 2D structure:', err);
    error.value = 'Failed to render 2D structure';
  }
};

const render3D = async () => {
  if (!canvas3d.value || !props.smiles) return;
  
  error.value = null;
  
  try {
    if (window.__rdkitReady) await window.__rdkitReady;
    if (!window.RDKit) {
      error.value = 'RDKit not available';
      return;
    }

    if (!window.$3Dmol) {
      error.value = '3Dmol.js not available';
      return;
    }

    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const RDKit = window.RDKit as any;
    const $3Dmol = window.$3Dmol;
    
    // Generate 3D coordinates using RDKit
    const mol = RDKit.get_mol(props.smiles);
    if (!mol || !mol.is_valid()) {
      error.value = 'Invalid SMILES structure';
      return;
    }

    // Add hydrogens and generate 3D conformer
    mol.add_hs();
    const molblock = mol.get_molblock();
    mol.delete();

    // Clear previous viewer if exists
    if (viewer3d) {
      viewer3d.clear();
    }

    // Create 3Dmol viewer
    canvas3d.value.innerHTML = '';
    viewer3d = $3Dmol.createViewer(canvas3d.value);
    
    // Set transparent background
    viewer3d.setBackgroundColor(0xffffff, 0);
    
    // Add molecule
    viewer3d.addModel(molblock, 'mol');
    viewer3d.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.3 } });
    viewer3d.zoomTo();
    viewer3d.render();
  } catch (err) {
    console.error('Error rendering 3D structure:', err);
    error.value = 'Failed to render 3D structure';
  }
};

const renderMolecule = async () => {
  await nextTick();
  
  if (props.mode === '2d') {
    await render2D();
  } else if (props.mode === '3d') {
    await render3D();
  }
};

watch(() => [props.smiles, props.mode], renderMolecule, { immediate: false });

onMounted(async () => {
  await renderMolecule();
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-molecule-renderer {
  width: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 1rem;

  &__canvas {
    width: 100%;
    height: auto;
    min-height: 200px;
    display: flex;
    align-items: center;
    justify-content: center;
    border-radius: 12px;
    background: color-with-opacity(var(--surface-white), $opacity-medium);
    overflow: hidden;
    
    &--3d {
      position: relative;
      aspect-ratio: 4 / 3;
    }
  }

  &__error {
    padding: 0.5rem 1rem;
    border-radius: 4px;
    background: rgba(255, 0, 0, 0.1);
  }
}

:deep(svg) {
  max-width: 100%;
  max-height: 100%;
}

:deep(svg rect) {
  fill: none !important;
}

:deep(svg path[fill="#000000"]) {
  fill: var(--text) !important;
}

:deep(svg path[style*="stroke:#000000"]) {
  stroke: var(--text) !important;
}

:deep(svg path[fill="#FF0000"]) {
  fill: var(--accent-red) !important;
}

:deep(svg path[style*="stroke:#FF0000"]) {
  stroke: var(--accent-red) !important;
}

:deep(svg path[fill="#0000FF"]) {
  fill: var(--accent-blue) !important;
}

:deep(svg path[style*="stroke:#0000FF"]) {
  stroke: var(--accent-blue) !important;
}
</style>
