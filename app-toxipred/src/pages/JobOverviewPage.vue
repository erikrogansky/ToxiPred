<template>
  <tp-page class="tp-job-overview column items-center text-center justify-center q-pa-md q-pa-lg-xl">

    <div class="tp-header full-width">
      <div class="tp-header__left">
        <tp-icon-button
          icon-name="arrow-left-02"
          weight="regular"
          :size="28"
          @click="$router.push('/workspace')"
        />
        <h1 class="tp-main-heading">{{ displayName }}</h1>
      </div>
      <div class="tp-header__right">
        <tp-button
          label="Share link"
          variant="outline"
          @click="showShareDialog = true"
        />
        <tp-button
          :label="generatingPdf ? 'Generating...' : 'Download PDF report'"
          :disabled="generatingPdf || !result"
          @click="downloadPdf"
        />
      </div>
    </div>

    <div v-if="loading">
      Loading prediction details…
    </div>

    <div v-else-if="error">
      {{ error }}
    </div>

    <div v-else-if="result" class="column items-center full-width">
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

      <div v-if="!smiles" class="text-caption q-mt-md">
        No structure available for visualization.
      </div>
    </div>

    <div v-else>
      No response available.
    </div>

    <!-- Share Dialog -->
    <tp-share-dialog v-model="showShareDialog" :job-id="jobId" />
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
import TpCompoundInfo from 'components/TpCompoundInfo.vue';
import TpShareDialog from 'components/TpShareDialog.vue';
import { api } from 'src/boot/axios';
import { useRoute } from 'vue-router';
import { ref, onMounted, computed } from 'vue';
import { useModelsStore } from 'src/stores/models-store';
import { jsPDF } from 'jspdf';

const route = useRoute();
const jobId = route.params.job_id?.toString() || '';
const modelsStore = useModelsStore();

// Dialog state
const showShareDialog = ref(false);
const generatingPdf = ref(false);

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

const loading = ref(true);
const error = ref<string | null>(null);
const result = ref<JobResultPayload | null>(null);
const smiles = ref<string>('');
const atomScores = ref<number[]>([]);
const featureNames = ref<string[]>([]);
const featureValues = ref<(number | null)[]>([]);
const featureScores = ref<number[]>([]);
const confidence = ref<number | null>(null);

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
    featureValues.value = result.value.feature_values ?? [];
    featureScores.value = (result.value.feature_scores ?? []);
    confidence.value = result.value.confidence ?? null;

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

const displayName = computed<string>(() => {
  if (!result.value) return 'Loading...';
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

// Descriptor type mapping for PDF
const descriptorTypes: Record<string, { type: string; displayName: string }> = {
  'qed': { type: 'Physicochemical', displayName: 'QED' },
  'BalabanJ': { type: 'Physicochemical', displayName: 'Balaban J' },
  'HallKierAlpha': { type: 'Physicochemical', displayName: 'Hall-Kier Alpha' },
  'MaxEStateIndex': { type: 'Electronic', displayName: 'Max EState Index' },
  'MinAbsEStateIndex': { type: 'Electronic', displayName: 'Min Abs EState Index' },
  'PEOE_VSA1': { type: 'Electronic', displayName: 'PEOE VSA1' },
  'PEOE_VSA3': { type: 'Electronic', displayName: 'PEOE VSA3' },
  'PEOE_VSA6': { type: 'Electronic', displayName: 'PEOE VSA6' },
  'PEOE_VSA7': { type: 'Electronic', displayName: 'PEOE VSA7' },
  'PEOE_VSA9': { type: 'Electronic', displayName: 'PEOE VSA9' },
  'PEOE_VSA11': { type: 'Electronic', displayName: 'PEOE VSA11' },
  'SMR_VSA6': { type: 'Electronic', displayName: 'SMR VSA6' },
  'SlogP_VSA2': { type: 'Electronic', displayName: 'SlogP VSA2' },
  'EState_VSA2': { type: 'Electronic', displayName: 'EState VSA2' },
  'EState_VSA3': { type: 'Electronic', displayName: 'EState VSA3' },
  'EState_VSA4': { type: 'Electronic', displayName: 'EState VSA4' },
  'VSA_EState2': { type: 'Electronic', displayName: 'VSA EState2' },
  'VSA_EState4': { type: 'Electronic', displayName: 'VSA EState4' },
  'VSA_EState6': { type: 'Electronic', displayName: 'VSA EState6' },
  'AvgIpc': { type: 'Topological', displayName: 'Avg IPC' },
  'BertzCT': { type: 'Topological', displayName: 'Bertz CT' },
  'Chi4n': { type: 'Topological', displayName: 'Chi4n' },
  'NumHeteroatoms': { type: 'Structural', displayName: 'Num Heteroatoms' },
  'fr_NH1': { type: 'Structural', displayName: 'NH1 Fragments' },
  'fr_NH2': { type: 'Structural', displayName: 'NH2 Fragments' },
  'fr_amide': { type: 'Structural', displayName: 'Amide Fragments' },
};

function formatDescriptorValue(value: number | null, name: string): string {
  if (value === null || Number.isNaN(value)) return '—';
  if (name.startsWith('Num') || name.startsWith('fr_')) {
    return Math.round(value).toString();
  }
  if (Math.abs(value) < 0.01 && value !== 0) {
    return value.toExponential(2);
  }
  if (Number.isInteger(value)) {
    return value.toString();
  }
  return value.toFixed(2);
}

// Type badge colors matching the app
const typeBadgeColors: Record<string, { bg: [number, number, number]; text: [number, number, number] }> = {
  'Physicochemical': { bg: [232, 180, 240], text: [92, 45, 107] },
  'Electronic': { bg: [145, 213, 245], text: [26, 95, 122] },
  'Structural': { bg: [168, 230, 207], text: [45, 92, 74] },
  'Topological': { bg: [255, 191, 135], text: [139, 69, 19] },
  'Other': { bg: [208, 208, 208], text: [74, 74, 74] },
};

// Helper to convert SVG to PNG data URL
async function svgToDataUrl(svgString: string, width: number, height: number): Promise<string | null> {
  return new Promise((resolve) => {
    const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(svgBlob);
    
    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement('canvas');
      canvas.width = width;
      canvas.height = height;
      const ctx = canvas.getContext('2d');
      if (ctx) {
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, width, height);
        ctx.drawImage(img, 0, 0, width, height);
        URL.revokeObjectURL(url);
        resolve(canvas.toDataURL('image/png'));
      } else {
        URL.revokeObjectURL(url);
        resolve(null);
      }
    };
    img.onerror = () => {
      URL.revokeObjectURL(url);
      resolve(null);
    };
    img.src = url;
  });
}
// Helper to render XSMILES visualization as image - uses exact same rendering as TpXsmilesRenderer
async function renderXsmilesImage(smilesStr: string, atomScores: number[], width: number, height: number): Promise<string | null> {
  try {
    if (window.__xsmilesReady) await window.__xsmilesReady;
    if (window.__rdkitReady) await window.__rdkitReady;
    if (!window.RDKit || !window.xsmiles) return null;

    // Create a temporary container (visible but off-screen for proper rendering)
    const tempDiv = document.createElement('div');
    tempDiv.style.cssText = `position: fixed; left: -9999px; top: 0; width: ${width}px; height: ${height}px; background: white;`;
    document.body.appendChild(tempDiv);

    // Use EXACT same setup as TpXsmilesRenderer component
    const setup = {
      molecule: {
        string: smilesStr,
        method: {
          name: 'atom_mask_dummy',
          scores: atomScores,
          attributes: {}
        },
        attributes: {}
      },
      width,
      height,
      bondLength: 75,
      gradientConfig: {
        palette: {
          name: 'MyPalette',
          colors: ['#7FFFD4', '#C7FFEA', '#FFFFFF', '#FCDDDD', '#F9B4B4']
        },
        colorDomain: [-0.09, -0.05, 0, 0.05, 0.09],
        highlight: true,
        blur: 0.7,
        opacity: { min: 0.7, max: 1.0 },
        radius: { min: 25, max: 50 },
        thresholds: [],
        delta: 0.05,
      },
      drawerType: 'RDKitDrawer_black',
      hideBarChart: true,
      hideAttributesTable: true,
      showScoresOnStructure: true
    };

    const api = window.xsmiles;
    const fn = api.appendSingleView || api.appendSingle || api.renderSingle || api.appendView;

    if (typeof fn !== 'function') {
      console.error('XSMILES API not found:', Object.keys(api));
      document.body.removeChild(tempDiv);
      return null;
    }

    try {
      fn(tempDiv, setup);
    } catch (e) {
      console.error('XSMILES render failed. Setup was:', setup, e);
      document.body.removeChild(tempDiv);
      return null;
    }

    // Wait for XSMILES to fully render (D3 animations + canvas drawing)
    await new Promise(resolve => setTimeout(resolve, 2500));

    // Find the SingleView container or use tempDiv
    const singleView = tempDiv.querySelector('.SingleView') as HTMLElement || tempDiv;
    
    // Manual screenshot approach: combine canvas and SVG elements
    // Find all canvases and SVGs in the container
    const canvases = Array.from(singleView.querySelectorAll('canvas'));
    const svgs = Array.from(singleView.querySelectorAll('svg'));
    
    if (canvases.length === 0 && svgs.length === 0) {
      console.error('No canvas or SVG elements found in XSMILES output');
      document.body.removeChild(tempDiv);
      return null;
    }

    // Calculate actual bounds of all elements
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    const containerRect = singleView.getBoundingClientRect();
    
    [...canvases, ...svgs].forEach(el => {
      const rect = el.getBoundingClientRect();
      const x = rect.left - containerRect.left;
      const y = rect.top - containerRect.top;
      minX = Math.min(minX, x);
      minY = Math.min(minY, y);
      maxX = Math.max(maxX, x + rect.width);
      maxY = Math.max(maxY, y + rect.height);
    });

    const actualWidth = maxX - minX;
    const actualHeight = maxY - minY;

    // Create a composite canvas with actual content size
    const compositeCanvas = document.createElement('canvas');
    compositeCanvas.width = actualWidth;
    compositeCanvas.height = actualHeight;
    const ctx = compositeCanvas.getContext('2d');
    
    if (!ctx) {
      document.body.removeChild(tempDiv);
      return null;
    }

    // Fill with white background
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, actualWidth, actualHeight);

    // Draw all canvas elements first (molecule structure)
    for (const canvas of canvases) {
      const rect = canvas.getBoundingClientRect();
      const x = rect.left - containerRect.left - minX;
      const y = rect.top - containerRect.top - minY;
      
      try {
        ctx.drawImage(canvas, x, y, rect.width, rect.height);
      } catch (err) {
        console.warn('Could not draw canvas (may be tainted):', err);
        // If canvas is tainted, we'll just skip it
      }
    }

    // Draw all SVG elements on top (text labels, scores) - WAIT for each to load
    for (const svg of svgs) {
      const rect = svg.getBoundingClientRect();
      const x = rect.left - containerRect.left - minX;
      const y = rect.top - containerRect.top - minY;
      
      // Convert SVG to image and draw it
      const svgString = new XMLSerializer().serializeToString(svg);
      const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
      const url = URL.createObjectURL(svgBlob);
      
      // Make sure we WAIT for the image to load before continuing
      const loaded = await new Promise<boolean>((resolve) => {
        const img = new Image();
        img.onload = () => {
          try {
            ctx.drawImage(img, x, y, rect.width, rect.height);
            console.log(`Drew SVG at (${x}, ${y}) with size ${rect.width}x${rect.height}`);
            resolve(true);
          } catch (err) {
            console.error('Error drawing SVG to canvas:', err);
            resolve(false);
          }
          URL.revokeObjectURL(url);
        };
        img.onerror = (err) => {
          console.error('Error loading SVG as image:', err);
          URL.revokeObjectURL(url);
          resolve(false);
        };
        img.src = url;
      });
      
      if (!loaded) {
        console.warn('Failed to draw SVG element');
      }
    }

    console.log(`Created composite canvas: ${actualWidth}x${actualHeight} (from bounds: minX=${minX}, minY=${minY}, maxX=${maxX}, maxY=${maxY})`);
    
    const imgData = compositeCanvas.toDataURL('image/png');
    document.body.removeChild(tempDiv);
    return imgData;
  } catch (err) {
    console.error('Error rendering XSMILES:', err);
    return null;
  }
}

// Helper to render 3D molecule as image
async function render3DMoleculeImage(smilesStr: string, width: number, height: number): Promise<string | null> {
  try {
    if (window.__rdkitReady) await window.__rdkitReady;
    if (!window.RDKit || !window.$3Dmol) return null;

    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const RDKit = window.RDKit as any;
    const $3Dmol = window.$3Dmol;

    // Generate molblock with 3D coordinates
    const mol = RDKit.get_mol(smilesStr);
    if (!mol || !mol.is_valid()) return null;

    mol.add_hs();
    const molblock = mol.get_molblock();
    mol.delete();

    // Create a temporary container
    const tempDiv = document.createElement('div');
    tempDiv.style.cssText = `position: fixed; left: -9999px; width: ${width}px; height: ${height}px; background: white;`;
    document.body.appendChild(tempDiv);

    // Create 3Dmol viewer
    const viewer = $3Dmol.createViewer(tempDiv, { backgroundColor: 'white' });
    viewer.addModel(molblock, 'mol');
    viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.3 } });
    viewer.zoomTo();
    viewer.render();

    // Wait for rendering to complete
    await new Promise(resolve => setTimeout(resolve, 300));

    // Get image from viewer
    const imgData = viewer.pngURI();
    
    // Cleanup
    viewer.clear();
    document.body.removeChild(tempDiv);

    return imgData;
  } catch (err) {
    console.error('Error rendering 3D molecule:', err);
    return null;
  }
}

// Helper to render waterfall chart as image
async function renderWaterfallChart(
  features: string[], 
  scores: number[], 
  width: number, 
  height: number
): Promise<string | null> {
  try {
    if (window.__plotlyReady) await window.__plotlyReady;
    const plotly = window.Plotly;
    if (!plotly) return null;
    
    // Prepare data
    const pairs: { name: string; score: number }[] = [];
    const len = Math.min(features.length, scores.length);
    
    for (let i = 0; i < len; i++) {
      const s = scores[i];
      if (s == null || Number.isNaN(s)) continue;
      pairs.push({ name: features[i] ?? '', score: s });
    }
    
    const sorted = pairs.sort((a, b) => Math.abs(b.score) - Math.abs(a.score)).slice(0, 15);
    if (!sorted.length) return null;
    
    const labels = sorted.map(f => f.name);
    const xValues = sorted.map(f => f.score);
    const colors = sorted.map(f => f.score >= 0 ? '#ef4444' : '#3b82f6');
    const text = sorted.map(f => f.score >= 0 ? `+${f.score.toFixed(3)}` : f.score.toFixed(3));
    
    // Create a temporary div
    const tempDiv = document.createElement('div');
    tempDiv.style.cssText = 'position: fixed; left: -9999px; width: ' + width + 'px; height: ' + height + 'px;';
    document.body.appendChild(tempDiv);
    
    const data: Partial<Plotly.PlotData>[] = [{
      type: 'bar',
      orientation: 'h',
      y: labels,
      x: xValues,
      text,
      textposition: 'outside',
      marker: { color: colors },
      hoverinfo: 'none'
    }];
    
    const layout: Partial<Plotly.Layout> = {
      margin: { l: 100, r: 60, t: 20, b: 40 },
      xaxis: {
        title: { text: 'Contribution to prediction' },
        zeroline: true,
        zerolinecolor: 'rgba(0,0,0,0.3)',
        zerolinewidth: 1
      },
      yaxis: { automargin: true },
      showlegend: false,
      paper_bgcolor: '#ffffff',
      plot_bgcolor: '#ffffff',
      font: { size: 10 }
    };
    
    await plotly.newPlot(tempDiv, data, layout, { displayModeBar: false, staticPlot: true });
    
    // Convert to image
    const imgData = await plotly.toImage(tempDiv, { format: 'png', width, height, scale: 2 });
    
    // Cleanup
    plotly.purge(tempDiv);
    document.body.removeChild(tempDiv);
    
    return imgData;
  } catch (err) {
    console.error('Error rendering waterfall chart:', err);
    return null;
  }
}

async function downloadPdf() {
  if (!result.value) return;
  
  generatingPdf.value = true;
  
  try {
    const doc = new jsPDF({
      orientation: 'portrait',
      unit: 'mm',
      format: 'a4'
    });
    
    const pageWidth = doc.internal.pageSize.getWidth();
    const pageHeight = doc.internal.pageSize.getHeight();
    const margin = 10;
    const contentWidth = pageWidth - 2 * margin;
    let y = margin;
    
    // Colors matching your app
    const primaryGreen: [number, number, number] = [23, 185, 136];      // #17B988
    const primaryGreenLight: [number, number, number] = [63, 219, 165]; // #3FDBA5
    const darkText: [number, number, number] = [41, 50, 48];            // #293230
    const grayText: [number, number, number] = [82, 98, 93];            // #52625D
    const lightGray: [number, number, number] = [232, 255, 243];        // #E8FFF3 (green-50)
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    const mediumGray: [number, number, number] = [218, 255, 237];       // #DAFFED
    const dangerRed: [number, number, number] = [239, 68, 68];          // #EF4444
    const confidenceYellow: [number, number, number] = [245, 200, 66];  // #F5C842
    const white: [number, number, number] = [255, 255, 255];
    
    // Helper to check and add new page
    const checkNewPage = (neededHeight: number) => {
      if (y + neededHeight > pageHeight - margin) {
        doc.addPage();
        y = margin;
        return true;
      }
      return false;
    };
    
    // Helper to draw rounded rectangle
    const drawRoundedRect = (x: number, yPos: number, w: number, h: number, r: number, fill: [number, number, number], stroke?: [number, number, number]) => {
      doc.setFillColor(...fill);
      if (stroke) {
        doc.setDrawColor(...stroke);
        doc.setLineWidth(0.3);
        doc.roundedRect(x, yPos, w, h, r, r, 'FD');
      } else {
        doc.roundedRect(x, yPos, w, h, r, r, 'F');
      }
    };
    
    // ===== HEADER =====
    // Helper to create logo as SVG and convert to image
    const createLogoImage = async (): Promise<string | null> => {
      const logoSvg = `
        <svg xmlns="http://www.w3.org/2000/svg" width="160" height="32" viewBox="0 0 160 32" fill="#17B988">
          <path d="M20.844 5.652H13.464V27.054H7.38V5.652H0V0.828H20.844V5.652Z"/>
          <path d="M27.6376 8.1C29.0536 8.1 30.3436 8.322 31.5076 8.766C32.6716 9.21 33.6676 9.846 34.4956 10.674C35.3356 11.502 35.9836 12.51 36.4396 13.698C36.9076 14.874 37.1416 16.2 37.1416 17.676C37.1416 19.164 36.9076 20.508 36.4396 21.708C35.9836 22.896 35.3356 23.91 34.4956 24.75C33.6676 25.578 32.6716 26.22 31.5076 26.676C30.3436 27.12 29.0536 27.342 27.6376 27.342C26.2096 27.342 24.9076 27.12 23.7316 26.676C22.5676 26.22 21.5596 25.578 20.7076 24.75C19.8676 23.91 19.2136 22.896 18.7456 21.708C18.2896 20.508 18.0616 19.164 18.0616 17.676C18.0616 16.2 18.2896 14.874 18.7456 13.698C19.2136 12.51 19.8676 11.502 20.7076 10.674C21.5596 9.846 22.5676 9.21 23.7316 8.766C24.9076 8.322 26.2096 8.1 27.6376 8.1ZM27.6376 23.202C28.9216 23.202 29.8636 22.752 30.4636 21.852C31.0756 20.94 31.3816 19.56 31.3816 17.712C31.3816 15.864 31.0756 14.49 30.4636 13.59C29.8636 12.69 28.9216 12.24 27.6376 12.24C26.3176 12.24 25.3516 12.69 24.7396 13.59C24.1276 14.49 23.8216 15.864 23.8216 17.712C23.8216 19.56 24.1276 20.94 24.7396 21.852C25.3516 22.752 26.3176 23.202 27.6376 23.202Z"/>
          <path d="M56.9287 27.054H51.5647C51.1687 27.054 50.8447 26.958 50.5927 26.766C50.3527 26.574 50.1607 26.346 50.0167 26.082L46.6867 20.358C46.6267 20.562 46.5547 20.754 46.4707 20.934C46.3867 21.114 46.2967 21.282 46.2007 21.438L43.4647 26.082C43.2967 26.322 43.0987 26.544 42.8707 26.748C42.6427 26.952 42.3487 27.054 41.9887 27.054H37.0207C37.0207 27.054 43.2847 18.651 43.2847 17.406C43.2847 16.161 42.589 14.8706 41.6953 13.3706C40.8016 11.8706 39.8019 10.2463 37.7888 8.77055C36.1623 7.5782 34.0279 6.53637 33.2681 6.18021C33.18 6.13888 33.1958 6.01306 33.2922 5.99931C34.0406 5.89251 36.0541 5.69704 38.2155 6.26647C40.9492 6.98667 42.8472 8.59374 44.4028 10.2394C45.9585 11.885 47.3347 14.58 47.3347 14.58C47.4067 14.364 47.4967 14.154 47.6047 13.95C47.7127 13.734 47.8327 13.518 47.9647 13.302L50.3407 9.144C50.5087 8.88 50.6827 8.688 50.8627 8.568C51.0547 8.448 51.2887 8.388 51.5647 8.388H56.6767C56.6767 8.388 50.6647 15.5415 50.6647 17.154C50.6647 18.7665 56.9287 27.054 56.9287 27.054Z"/>
          <path d="M64.3607 8.388V27.054H58.7807V8.388H64.3607ZM64.9367 3.276C64.9367 3.72 64.8467 4.14 64.6667 4.536C64.4867 4.92 64.2407 5.262 63.9287 5.562C63.6167 5.85 63.2507 6.084 62.8307 6.264C62.4227 6.432 61.9847 6.516 61.5167 6.516C61.0607 6.516 60.6347 6.432 60.2387 6.264C59.8427 6.084 59.4887 5.85 59.1767 5.562C58.8767 5.262 58.6367 4.92 58.4567 4.536C58.2887 4.14 58.2047 3.72 58.2047 3.276C58.2047 2.82 58.2887 2.394 58.4567 1.998C58.6367 1.602 58.8767 1.254 59.1767 0.954C59.4887 0.654 59.8427 0.42 60.2387 0.251999C60.6347 0.0839996 61.0607 0 61.5167 0C61.9847 0 62.4227 0.0839996 62.8307 0.251999C63.2507 0.42 63.6167 0.654 63.9287 0.954C64.2407 1.254 64.4867 1.602 64.6667 1.998C64.8467 2.394 64.9367 2.82 64.9367 3.276Z"/>
          <path d="M77.9062 13.698C79.4902 13.698 80.6302 13.32 81.3262 12.564C82.0222 11.808 82.3703 10.752 82.3703 9.396C82.3703 8.796 82.2803 8.25 82.1003 7.758C81.9202 7.266 81.6443 6.846 81.2723 6.498C80.9123 6.138 80.4503 5.862 79.8863 5.67C79.3343 5.478 78.6742 5.382 77.9062 5.382H74.7383V13.698H77.9062ZM77.9062 0.828C79.7542 0.828 81.3382 1.05 82.6582 1.494C83.9902 1.926 85.0823 2.526 85.9342 3.294C86.7982 4.062 87.4342 4.968 87.8422 6.012C88.2502 7.056 88.4543 8.184 88.4543 9.396C88.4543 10.704 88.2442 11.904 87.8242 12.996C87.4043 14.088 86.7623 15.024 85.8983 15.804C85.0343 16.584 83.9362 17.196 82.6042 17.64C81.2842 18.072 79.7182 18.288 77.9062 18.288H74.7383V27.054H68.6543V0.828H77.9062Z"/>
          <path d="M95.8669 11.268C96.4429 10.26 97.1029 9.468 97.8469 8.892C98.6029 8.316 99.4669 8.028 100.439 8.028C101.279 8.028 101.963 8.226 102.491 8.622L102.131 12.726C102.071 12.99 101.969 13.17 101.825 13.266C101.693 13.362 101.507 13.41 101.267 13.41C101.171 13.41 101.045 13.404 100.889 13.392C100.733 13.38 100.571 13.368 100.403 13.356C100.235 13.332 100.061 13.314 99.8809 13.302C99.7129 13.278 99.5569 13.266 99.4129 13.266C98.5969 13.266 97.9429 13.482 97.4509 13.914C96.9709 14.346 96.5629 14.946 96.2269 15.714V27.054H90.6469V8.388H93.9589C94.2349 8.388 94.4629 8.412 94.6429 8.46C94.8349 8.508 94.9909 8.586 95.1109 8.694C95.2429 8.79 95.3389 8.922 95.3989 9.09C95.4709 9.258 95.5309 9.462 95.5789 9.702L95.8669 11.268Z"/>
          <path d="M115.687 15.48C115.687 15.036 115.627 14.604 115.507 14.184C115.399 13.764 115.213 13.392 114.949 13.068C114.697 12.732 114.361 12.462 113.941 12.258C113.521 12.054 113.011 11.952 112.411 11.952C111.355 11.952 110.527 12.258 109.927 12.87C109.327 13.47 108.931 14.34 108.739 15.48H115.687ZM108.667 18.72C108.847 20.28 109.321 21.414 110.089 22.122C110.857 22.818 111.847 23.166 113.059 23.166C113.707 23.166 114.265 23.088 114.733 22.932C115.201 22.776 115.615 22.602 115.975 22.41C116.347 22.218 116.683 22.044 116.983 21.888C117.295 21.732 117.619 21.654 117.955 21.654C118.399 21.654 118.735 21.816 118.963 22.14L120.583 24.138C120.007 24.798 119.383 25.338 118.711 25.758C118.051 26.166 117.367 26.49 116.659 26.73C115.963 26.958 115.261 27.114 114.553 27.198C113.857 27.294 113.191 27.342 112.555 27.342C111.247 27.342 110.017 27.132 108.865 26.712C107.725 26.28 106.723 25.644 105.859 24.804C105.007 23.964 104.329 22.92 103.825 21.672C103.333 20.424 103.087 18.972 103.087 17.316C103.087 16.068 103.297 14.886 103.717 13.77C104.149 12.654 104.761 11.676 105.553 10.836C106.357 9.996 107.323 9.33 108.451 8.838C109.591 8.346 110.875 8.1 112.303 8.1C113.539 8.1 114.667 8.292 115.687 8.676C116.719 9.06 117.601 9.618 118.333 10.35C119.077 11.082 119.653 11.982 120.061 13.05C120.481 14.106 120.691 15.306 120.691 16.65C120.691 17.07 120.673 17.412 120.637 17.676C120.601 17.94 120.535 18.15 120.439 18.306C120.343 18.462 120.211 18.57 120.043 18.63C119.887 18.69 119.683 18.72 119.431 18.72H108.667Z"/>
          <path d="M134.805 13.392C134.385 12.924 133.929 12.6 133.437 12.42C132.945 12.228 132.429 12.132 131.889 12.132C131.373 12.132 130.899 12.234 130.467 12.438C130.035 12.642 129.657 12.972 129.333 13.428C129.021 13.872 128.775 14.454 128.595 15.174C128.415 15.894 128.325 16.764 128.325 17.784C128.325 18.768 128.397 19.596 128.541 20.268C128.685 20.94 128.883 21.486 129.135 21.906C129.399 22.314 129.711 22.608 130.071 22.788C130.431 22.968 130.833 23.058 131.277 23.058C131.685 23.058 132.051 23.022 132.375 22.95C132.711 22.866 133.017 22.752 133.293 22.608C133.569 22.464 133.827 22.284 134.067 22.068C134.319 21.84 134.565 21.582 134.805 21.294V13.392ZM140.385 0.107999V27.054H136.929C136.233 27.054 135.777 26.742 135.561 26.118L135.129 24.696C134.757 25.092 134.367 25.452 133.959 25.776C133.551 26.1 133.107 26.382 132.627 26.622C132.159 26.85 131.649 27.024 131.097 27.144C130.557 27.276 129.969 27.342 129.333 27.342C128.361 27.342 127.461 27.126 126.633 26.694C125.805 26.262 125.091 25.638 124.491 24.822C123.891 24.006 123.417 23.01 123.069 21.834C122.733 20.646 122.565 19.296 122.565 17.784C122.565 16.392 122.757 15.102 123.141 13.914C123.525 12.714 124.065 11.682 124.761 10.818C125.469 9.942 126.315 9.258 127.299 8.766C128.283 8.274 129.369 8.028 130.557 8.028C131.529 8.028 132.345 8.166 133.005 8.442C133.665 8.718 134.265 9.09 134.805 9.558V0.107999H140.385Z"/>
        </svg>
      `;
      return await svgToDataUrl(logoSvg, 160, 32);
    };
    
    // Logo
    const logoImg = await createLogoImage();
    if (logoImg) {
      doc.addImage(logoImg, 'PNG', margin, y, 38, 7.6);
      //y += 7.6;
    } else {
      drawRoundedRect(margin, y, 38, 10, 2, lightGray);
      doc.setDrawColor(...primaryGreen);
      doc.setLineWidth(0.6);
      doc.circle(margin + 6, y + 5, 3.5);
      doc.setLineWidth(0.5);
      doc.line(margin + 4, y + 5, margin + 5.5, y + 6.5);
      doc.line(margin + 5.5, y + 6.5, margin + 8, y + 3.5);
      doc.setFontSize(12);
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...darkText);
      doc.text('ToxiPred', margin + 12, y + 6.5);
      //y += 10;
    }
    
    // Meta info (right side)
    doc.setFontSize(8);
    doc.setFont('helvetica', 'normal');
    doc.setTextColor(...grayText);
    const currentDate = new Date().toLocaleDateString('en-US', {
      year: 'numeric', month: 'long', day: 'numeric',
      hour: '2-digit', minute: '2-digit'
    });
    doc.text(`Report ID: ${jobId}`, pageWidth - margin, y + 3, { align: 'right' });
    doc.text(`Generated: ${currentDate}`, pageWidth - margin, y + 7, { align: 'right' });
    
    y += 24;
    
    // ===== TITLE =====
    doc.setFontSize(20);
    doc.setFont('helvetica', 'bold');
    doc.setTextColor(...darkText);
    doc.text('Toxicity Prediction Report', pageWidth / 2, y, { align: 'center' });
    y += 12;
    
    // ===== PREDICTION RESULT BADGES (matching app style) =====
    const badgeHeight = 12;
    const badgeGap = 6;
    const badgeWidth = (contentWidth - badgeGap) / 2;
    
    const isToxic = predictionValue.value === 1;
    const resultBgColor = isToxic ? dangerRed : primaryGreenLight;
    const resultTextColor = isToxic ? white : darkText;
    
    // Result badge
    drawRoundedRect(margin, y, badgeWidth, badgeHeight, 2, resultBgColor);
    doc.setFontSize(11);
    doc.setFont('helvetica', 'bold');
    doc.setTextColor(...resultTextColor);
    const resultLabel = predictionValue.value === 0 ? negativeLabel.value : positiveLabel.value;
    doc.text(resultLabel, margin + badgeWidth / 2, y + 7.5, { align: 'center' });
    
    // Confidence badge
    drawRoundedRect(margin + badgeWidth + badgeGap, y, badgeWidth, badgeHeight, 2, confidenceYellow);
    doc.setTextColor(...darkText);
    const confidenceText = confidence.value !== null ? `Confidence: ${(confidence.value * 100).toFixed(1)}%` : 'Confidence: —';
    doc.text(confidenceText, margin + badgeWidth + badgeGap + badgeWidth / 2, y + 7.5, { align: 'center' });
    
    y += badgeHeight + 10;
    
    // ===== COMPOUND INFORMATION =====
    const compoundInfo: [string, string][] = [];
    if (result.value.trivial_name) compoundInfo.push(['Trivial Name', result.value.trivial_name]);
    if (result.value.name && result.value.name !== result.value.trivial_name) {
      compoundInfo.push(['Chemical Name', result.value.name]);
    }
    if (result.value.formula) compoundInfo.push(['Molecular Formula', result.value.formula]);
    if (smiles.value) compoundInfo.push(['SMILES', smiles.value]);
    if (result.value.input_type) {
      const typeMap: Record<string, string> = {
        'smiles': 'SMILES', 'inchi': 'InChI', 'cas': 'CAS Number', 'name': 'Chemical Name'
      };
      compoundInfo.push(['Input Type', typeMap[result.value.input_type] || result.value.input_type]);
    }
    if (result.value.other_names && result.value.other_names.length > 0) {
      compoundInfo.push(['Other Names', result.value.other_names.slice(0, 3).join(', ')]);
    }
    
    // Calculate content height
    let contentHeight = 0;
    for (const [, value] of compoundInfo) {
      const lines = doc.splitTextToSize(value, contentWidth - 50);
      contentHeight += lines.length * 4.5 + 2;
    }
    
    checkNewPage(contentHeight + 20);
    
    // Section header
    doc.setFontSize(12);
    doc.setFont('helvetica', 'bold');
    doc.setTextColor(...darkText);
    doc.text('Compound Information', margin, y);
    y += 2;
    doc.setDrawColor(...primaryGreen);
    doc.setLineWidth(0.5);
    doc.line(margin, y, margin + 42, y);
    y += 8;
    
    // Compound info rows
    doc.setFontSize(9);
    for (const [label, value] of compoundInfo) {
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...grayText);
      doc.text(label + ':', margin, y);
      
      doc.setFont('helvetica', 'normal');
      doc.setTextColor(...darkText);
      
      const maxValueWidth = contentWidth - 50;
      const splitValue = doc.splitTextToSize(value, maxValueWidth);
      doc.text(splitValue, margin + 45, y);
      y += splitValue.length * 4.5 + 2;
    }
    
    y += 8;
    
    // ===== MOLECULAR STRUCTURE (2D and 3D side by side) =====
    if (smiles.value) {
      const structureWidth = (contentWidth - 6) / 2; // Split width with gap
      const structureHeight = structureWidth; // Make containers square for proper aspect ratio
      
      checkNewPage(structureHeight + 20);
      
      // Section header
      doc.setFontSize(12);
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...darkText);
      doc.text('Molecular Structure', margin, y);
      y += 2;
      doc.setDrawColor(...primaryGreen);
      doc.setLineWidth(0.5);
      doc.line(margin, y, margin + 40, y);
      y += 6;
      
      const structureY = y;
      
      // Render 2D structure (left side)
      try {
        if (window.__rdkitReady) await window.__rdkitReady;
        if (window.RDKit) {
          // eslint-disable-next-line @typescript-eslint/no-explicit-any
          const RDKit = window.RDKit as any;
          const mol = RDKit.get_mol(smiles.value);
          
          if (mol && mol.is_valid()) {
            const svgString = mol.get_svg(500, 300);
            mol.delete();
            
            const imgData = await svgToDataUrl(svgString, 500, 300);
            if (imgData) {
              // 2D label
              doc.setFontSize(9);
              doc.setFont('helvetica', 'bold');
              doc.setTextColor(...grayText);
              doc.text('2D Structure', margin + structureWidth / 2, structureY + 3, { align: 'center' });
              
              // 2D image - fit within container maintaining aspect ratio
              const img2dAspectRatio = 500 / 300; // width/height of original
              let img2dWidth = structureWidth;
              let img2dHeight = img2dWidth / img2dAspectRatio;
              
              // If height exceeds container, scale down
              if (img2dHeight > structureHeight) {
                img2dHeight = structureHeight;
                img2dWidth = img2dHeight * img2dAspectRatio;
              }
              
              const img2dX = margin + (structureWidth - img2dWidth) / 2; // Center horizontally
              const img2dY = structureY + 7 + (structureHeight - img2dHeight) / 2; // Center vertically
              
              drawRoundedRect(margin - 2, structureY + 5, structureWidth + 4, structureHeight + 4, 3, [250, 250, 250], [220, 220, 220]);
              doc.addImage(imgData, 'PNG', img2dX, img2dY, img2dWidth, img2dHeight);
            }
          }
        }
      } catch (err) {
        console.error('Error rendering 2D structure for PDF:', err);
      }
      
      // Render 3D structure (right side)
      const mol3dImg = await render3DMoleculeImage(smiles.value, 500, 500);
      if (mol3dImg) {
        // 3D label
        doc.setFontSize(9);
        doc.setFont('helvetica', 'bold');
        doc.setTextColor(...grayText);
        doc.text('3D Structure', margin + structureWidth + 6 + structureWidth / 2, structureY + 3, { align: 'center' });
        
        // 3D image - fit within container maintaining aspect ratio (already square)
        const img3dWidth = structureWidth;
        const img3dHeight = structureWidth;
        const img3dX = margin + structureWidth + 6;
        const img3dY = structureY + 7;
        
        drawRoundedRect(margin + structureWidth + 4, structureY + 5, structureWidth + 4, structureHeight + 4, 3, [250, 250, 250], [220, 220, 220]);
        doc.addImage(mol3dImg, 'PNG', img3dX, img3dY, img3dWidth, img3dHeight);
      }
      
      y += structureHeight + 16;
    }
    
    // ===== ATOMIC CONTRIBUTIONS (XSMILES-like) =====
    if (smiles.value && atomScores.value.length > 0) {
      checkNewPage(70);
      
      // Section header
      doc.setFontSize(12);
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...darkText);
      doc.text('Atomic Contributions', margin, y);
      y += 2;
      doc.setDrawColor(...primaryGreen);
      doc.setLineWidth(0.5);
      doc.line(margin, y, margin + 38, y);
      y += 4;
      
      // Description
      doc.setFontSize(8);
      doc.setFont('helvetica', 'normal');
      doc.setTextColor(...grayText);
      doc.text('Atoms colored by contribution: red = increases toxicity, green = decreases toxicity', margin, y);
      y += 6;
      
      // Use XSMILES renderer (same as UI component)
      const xsmilesImg = await renderXsmilesImage(smiles.value, atomScores.value, 800, 400);
      if (xsmilesImg) {
        const imgWidth = contentWidth;
        const imgHeight = contentWidth / 2; // Maintain 2:1 aspect ratio
        const imgX = margin;
        
        drawRoundedRect(imgX - 2, y - 2, imgWidth + 4, imgHeight + 4, 3, [250, 250, 250], [220, 220, 220]);
        doc.addImage(xsmilesImg, 'PNG', imgX, y, imgWidth, imgHeight);
        y += imgHeight + 8;
      }
    }
    
    // ===== DESCRIPTOR CONTRIBUTIONS CHART =====
    if (featureNames.value.length > 0 && featureScores.value.length > 0) {
      checkNewPage(85);
      
      // Section header
      doc.setFontSize(12);
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...darkText);
      doc.text('Descriptor Contributions', margin, y);
      y += 2;
      doc.setDrawColor(...primaryGreen);
      doc.setLineWidth(0.5);
      doc.line(margin, y, margin + 45, y);
      y += 4;
      
      // Description
      doc.setFontSize(8);
      doc.setFont('helvetica', 'normal');
      doc.setTextColor(...grayText);
      doc.text('Top descriptors contributing to the prediction (red = toxic, blue = non-toxic contribution)', margin, y);
      y += 6;
      
      // Render waterfall chart
      const chartImg = await renderWaterfallChart(featureNames.value, featureScores.value, 600, 350);
      if (chartImg) {
        // Maintain aspect ratio (600:350 from original) but fill width
        const chartAspectRatio = 600 / 350;
        const imgWidth = contentWidth;
        const imgHeight = imgWidth / chartAspectRatio;
        
        drawRoundedRect(margin - 2, y - 2, imgWidth + 4, imgHeight + 4, 3, white, [220, 220, 220]);
        doc.addImage(chartImg, 'PNG', margin, y, imgWidth, imgHeight);
        y += imgHeight + 10;
      }
    }
    
    // ===== MOLECULAR DESCRIPTORS TABLE =====
    if (featureNames.value.length > 0) {
      checkNewPage(50);
      
      // Section header
      doc.setFontSize(12);
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...darkText);
      doc.text('Molecular Descriptors', margin, y);
      y += 2;
      doc.setDrawColor(...primaryGreen);
      doc.setLineWidth(0.5);
      doc.line(margin, y, margin + 42, y);
      y += 6;
      
      // Table header
      const colWidths = { desc: 50, type: 32, value: 28, importance: 35 };
      const tableX = margin + 2;
      
      drawRoundedRect(margin, y, contentWidth, 8, 2, [235, 235, 235]);
      
      doc.setFontSize(8);
      doc.setFont('helvetica', 'bold');
      doc.setTextColor(...darkText);
      doc.text('Descriptor', tableX + 2, y + 5.5);
      doc.text('Type', tableX + colWidths.desc + 2, y + 5.5);
      doc.text('Value', tableX + colWidths.desc + colWidths.type + 2, y + 5.5);
      doc.text('Importance', tableX + colWidths.desc + colWidths.type + colWidths.value + 2, y + 5.5);
      
      y += 10;
      
      // Calculate max score for star rating
      let maxAbsScore = 0;
      for (const score of featureScores.value) {
        if (score != null && !Number.isNaN(score)) {
          maxAbsScore = Math.max(maxAbsScore, Math.abs(score));
        }
      }
      
      // Sort and build descriptor data
      const descriptorData: { displayName: string; type: string; value: string; stars: number }[] = [];
      const len = Math.min(featureNames.value.length, featureValues.value.length, featureScores.value.length);
      
      for (let i = 0; i < len; i++) {
        const name = featureNames.value[i] ?? '';
        const value = featureValues.value[i] ?? null;
        const score = featureScores.value[i] ?? 0;
        
        const info = descriptorTypes[name] || { type: 'Other', displayName: name };
        const normalizedScore = maxAbsScore > 0 ? Math.abs(score) / maxAbsScore : 0;
        const stars = Math.max(1, Math.min(5, Math.ceil(normalizedScore * 5)));
        
        descriptorData.push({
          displayName: info.displayName,
          type: info.type,
          value: formatDescriptorValue(value, name),
          stars
        });
      }
      
      descriptorData.sort((a, b) => b.stars - a.stars);
      
      // Table rows
      doc.setFont('helvetica', 'normal');
      let rowIndex = 0;
      for (const desc of descriptorData) {
        if (checkNewPage(7)) {
          // Redraw header on new page
          drawRoundedRect(margin, y, contentWidth, 8, 2, [235, 235, 235]);
          doc.setFontSize(10);
          doc.setFont('helvetica', 'bold');
          doc.setTextColor(...darkText);
          doc.text('Descriptor', tableX + 2, y + 5.5);
          doc.text('Type', tableX + colWidths.desc + 2, y + 5.5);
          doc.text('Value', tableX + colWidths.desc + colWidths.type + 2, y + 5.5);
          doc.text('Importance', tableX + colWidths.desc + colWidths.type + colWidths.value + 2, y + 5.5);
          y += 10;
          doc.setFont('helvetica', 'normal');
        }
        
        // Alternate row background
        if (rowIndex % 2 === 0) {
          doc.setFillColor(248, 250, 249);
          doc.rect(margin, y - 1, contentWidth, 6.5, 'F');
        }
        
        doc.setFontSize(10);
        doc.setTextColor(...darkText);
        doc.text(desc.displayName, tableX + 2, y + 3.5);
        
        // Type badge
        const badgeColors = typeBadgeColors[desc.type] ?? typeBadgeColors['Other']!;
        const typeText = desc.type;
        doc.setFontSize(8);
        doc.setFont('helvetica', 'bold');
        const typeBadgeWidth = doc.getTextWidth(typeText) + 5;
        drawRoundedRect(tableX + colWidths.desc, y - 0.5, typeBadgeWidth, 5, 1.5, badgeColors.bg);
        doc.setTextColor(...badgeColors.text);
        doc.text(typeText, tableX + colWidths.desc + 2.5, y + 3);
        
        // Value
        doc.setFontSize(10);
        doc.setFont('helvetica', 'normal');
        doc.setTextColor(...darkText);
        doc.text(desc.value, tableX + colWidths.desc + colWidths.type + 2, y + 3.5);
        
        // Star rating - draw SVG stars
        const starX = tableX + colWidths.desc + colWidths.type + colWidths.value + 2;
        const starY = y - 0.2;
        const starSize = 5;
        const starGap = 0.5;
        
        for (let s = 0; s < 5; s++) {
          const x = starX + s * (starSize + starGap);
          const isFilled = s < desc.stars;
          
          if (isFilled) {
            doc.setFillColor(...primaryGreenLight);
          } else {
            doc.setFillColor(200, 200, 200);
          }
          
          // Draw 5-pointed star
          const cx = x + starSize / 2;
          const cy = starY + starSize / 2;
          const outerRadius = starSize / 2;
          const innerRadius = starSize / 5;
          
          const points: [number, number][] = [];
          for (let i = 0; i < 5; i++) {
            const outerAngle = (i * 2 * Math.PI / 5) - Math.PI / 2;
            points.push([cx + outerRadius * Math.cos(outerAngle), cy + outerRadius * Math.sin(outerAngle)]);
            const innerAngle = ((i + 0.5) * 2 * Math.PI / 5) - Math.PI / 2;
            points.push([cx + innerRadius * Math.cos(innerAngle), cy + innerRadius * Math.sin(innerAngle)]);
          }
          
          doc.setDrawColor(0, 0, 0, 0); // No stroke
          // eslint-disable-next-line @typescript-eslint/no-explicit-any
          (doc as any).lines(
            points.slice(1).map((p, i) => [p[0] - points[i]![0], p[1] - points[i]![1]]),
            points[0]![0],
            points[0]![1],
            [1, 1],
            'F'
          );
        }
        
        doc.setFont('helvetica', 'normal');
        
        y += 6.5;
        rowIndex++;
      }
      
      y += 8;
    }
    
    // ===== FOOTER =====
    checkNewPage(20);
    
    doc.setDrawColor(...primaryGreen);
    doc.setLineWidth(0.5);
    doc.line(margin, y, pageWidth - margin, y);
    y += 6;
    
    doc.setFontSize(8);
    doc.setFont('helvetica', 'normal');
    doc.setTextColor(...grayText);
    doc.text(`Generated by ToxiPred • ${currentDate}`, pageWidth / 2, y, { align: 'center' });
    y += 5;
    
    doc.setFontSize(7);
    doc.setTextColor(150, 150, 150);
    const disclaimer = 'This report is for informational purposes only. Predictions are based on computational models and should be validated experimentally.';
    const disclaimerLines = doc.splitTextToSize(disclaimer, contentWidth);
    doc.text(disclaimerLines, pageWidth / 2, y, { align: 'center' });
    
    // Save PDF
    const filename = `toxipred-report-${displayName.value.replace(/[^a-zA-Z0-9]/g, '-').toLowerCase()}.pdf`;
    doc.save(filename);
    
  } catch (err) {
    console.error('Error generating PDF:', err);
  } finally {
    generatingPdf.value = false;
  }
}
</script>

<style scoped lang="scss">
.tp-job-overview {
  --content-gap: 32px;
  
  @media (min-width: 1024px) {
    --content-gap: 48px;
  }
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
</style>
