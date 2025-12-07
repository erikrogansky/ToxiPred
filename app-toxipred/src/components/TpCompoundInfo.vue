<template>
  <div class="tp-compound-info">
    <div class="tp-status-banner row items-center justify-center">
      <tp-icon 
        iconName="tick-circle" 
        :size="22"
        weight="regular"
      />
      <span class="paragraph black">{{ trivialName ? 'Matches known compound' : 'Unknown compound' }}</span>
    </div>

    <div class="tp-info-row">
      <div class="tp-info-label">Trivial name</div>
      <div class="tp-info-value">{{ displayName }}</div>
    </div>

    <div v-if="formula" class="tp-info-row">
      <div class="tp-info-label">Formula</div>
      <div class="tp-info-value" v-html="formattedFormula"></div>
    </div>

    <div v-if="otherNames && otherNames.length > 0" class="tp-info-row">
      <div class="tp-info-label">Other names</div>
      <div class="tp-info-value">{{ otherNames.slice(0, 5).join(', ') }}</div>
    </div>

    <div v-if="molecularWeight" class="tp-info-row">
      <div class="tp-info-label">Molecular weight</div>
      <div class="tp-info-value">{{ molecularWeight }} g/mol</div>
    </div>

    <div v-if="smiles" class="tp-info-row">
      <div class="tp-info-label">SMILES</div>
      <div class="tp-info-value tp-smiles">{{ smiles }}</div>
    </div>

    <div v-if="inputType" class="tp-info-row">
      <div class="tp-info-label">Input type</div>
      <div class="tp-info-value">{{ formattedInputType }}</div>
    </div>
  </div>
</template>

<script setup lang="ts">
import TpIcon from './TpIcon.vue';

import { computed } from 'vue';

interface Props {
  trivialName?: string | null | undefined;
  name?: string | null | undefined;
  formula?: string | null | undefined;
  otherNames?: string[] | null | undefined;
  smiles?: string | null | undefined;
  inputType?: string | null | undefined;
}

const props = defineProps<Props>();

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

const displayName = computed(() => {
  if (props.trivialName) return props.trivialName;
  if (props.smiles) return props.smiles;
  return props.name || 'Unknown compound';
});

const formattedFormula = computed(() => {
  if (!props.formula) return '';
  return props.formula.replace(/\d/g, digit => subscriptMap[digit] ?? digit);
});

const formattedInputType = computed(() => {
  if (!props.inputType) return '';
  const typeMap: Record<string, string> = {
    'smiles': 'SMILES',
    'inchi': 'InChI',
    'cas': 'CAS Number',
    'name': 'Chemical Name'
  };
  return typeMap[props.inputType] || props.inputType;
});

const molecularWeight = computed(() => {
  if (!props.formula) return null;
  
  try {
    const atomWeights: Record<string, number> = {
      'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
      'F': 18.998, 'P': 30.974, 'S': 32.06, 'Cl': 35.45,
      'Br': 79.904, 'I': 126.904, 'B': 10.81, 'Si': 28.085
    };
    
    let weight = 0;
    const matches = props.formula.matchAll(/([A-Z][a-z]?)(\d*)/g);
    
    for (const match of matches) {
      const element = match[1];
      const count = match[2] ? parseInt(match[2]) : 1;
      if (element && atomWeights[element]) {
        weight += atomWeights[element] * count;
      }
    }
    
    return weight > 0 ? weight.toFixed(2) : null;
  } catch {
    return null;
  }
});
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins.scss' as *;

.tp-compound-info {
  background: color-with-opacity(var(--surface-white), $opacity-medium);
  border-radius: 12px;
  overflow: hidden;
  max-width: 800px;
  width: 100%;
  padding: 14px 18px;
}

.tp-status-banner {
  display: flex;
  align-items: center;
  gap: 10px;
  padding: 10px 20px;
  background-color: var(--surface-brand-medium);
  font-weight: 500;
  font-size: 15px;
  border-radius: 6px;
}

.tp-info-row {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 16px 20px;
  border-bottom: 1px solid #E0E0E0;
  gap: 20px;
  
  &:last-child {
    border-bottom: none;
  }
}

.tp-info-label {
  font-weight: 600;
  color: #424242;
  font-size: 14px;
  flex-shrink: 0;
}

.tp-info-value {
  color: #212121;
  text-align: right;
  word-break: break-word;
  font-size: 14px;
  
  &.tp-smiles {
    font-family: 'Courier New', monospace;
    font-size: 13px;
  }
}
</style>
