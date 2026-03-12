<template>
  <tp-dialog
    :model-value="modelValue"
    title="Share Prediction"
    :loading="loading"
    loading-text="Creating share link..."
    :error="error"
    @update:model-value="emit('update:modelValue', $event)"
    @close="close"
    @retry="resetToConfig"
  >
    <!-- Initial State: Configure Password -->
    <div v-if="!shareData" class="tp-config-section">
      <p class="tp-description">
        Create a secure link to share this prediction. Choose to generate a random password or set your own.
      </p>

      <!-- Password Mode Toggle -->
      <tp-button-group
        :labels="['Generate Password', 'Custom Password']"
        :active-index="passwordMode === 'generate' ? 0 : 1"
        @click="passwordMode = $event === 0 ? 'generate' : 'custom'"
        class="q-mb-lg"
      />

      <!-- Custom Password Input -->
      <div v-if="passwordMode === 'custom'" class="tp-custom-password">
        <label class="tp-label">Your Password</label>
        <tp-glass-input
          v-model="customPassword"
          type="password"
          placeholder="Enter a secure password..."
        />
        <span class="tp-hint">
          <tp-icon icon-name="info-circle" weight="regular" :size="12" />
          Minimum 4 characters recommended
        </span>
      </div>

      <div v-else class="tp-generate-info">
        <div class="tp-info-glass">
          <tp-icon icon-name="shield-tick" weight="regular" :size="20" />
          <span>A secure 16-character password will be generated automatically</span>
        </div>
      </div>

      <div class="tp-actions">
        <tp-button
          label="Cancel"
          variant="outline"
          @click="close"
        />
        <tp-button
          label="Create Share Link"
          :disabled="passwordMode === 'custom' && customPassword.length < 1"
          @click="generateShareLink"
        />
      </div>
    </div>

    <!-- Success State -->
    <div v-else class="tp-share-success">
      <div class="tp-success-badge">
        <div class="tp-badge-icon">
          <tp-icon icon-name="tick-circle" weight="regular" :size="24" />
        </div>
        <span>Link Created Successfully</span>
      </div>

      <div class="tp-share-field">
        <label class="tp-label">
          <tp-icon icon-name="link" weight="regular" :size="14" />
          Share Link
        </label>
        <div class="tp-copy-field">
          <tp-glass-input
            :model-value="shareUrl"
            readonly
          />
          <button
            class="tp-copy-btn"
            :class="{ copied: copiedLink }"
            @click="copyLink"
          >
            <tp-icon :icon-name="copiedLink ? 'tick-circle' : 'copy'" weight="regular" :size="18" />
            <span>{{ copiedLink ? 'Copied!' : 'Copy' }}</span>
          </button>
        </div>
      </div>

      <div class="tp-share-field">
        <label class="tp-label">
          <tp-icon icon-name="key" weight="regular" :size="14" />
          Password
        </label>
        <div class="tp-copy-field">
          <tp-glass-input
            :model-value="shareData.password"
            readonly
            monospace
          />
          <button
            class="tp-copy-btn"
            :class="{ copied: copiedPassword }"
            @click="copyPassword"
          >
            <tp-icon :icon-name="copiedPassword ? 'tick-circle' : 'copy'" weight="regular" :size="18" />
            <span>{{ copiedPassword ? 'Copied!' : 'Copy' }}</span>
          </button>
        </div>
        <span class="tp-hint">
          <tp-icon icon-name="shield-tick" weight="regular" :size="12" />
          Share this password separately with the recipient
        </span>
      </div>

      <div class="tp-expiry-notice">
        <tp-icon icon-name="clock" weight="regular" :size="16" />
        <span>Expires {{ formatExpiry(shareData.expires_at) }}</span>
      </div>

      <div class="tp-actions">
        <tp-button
          label="Create Another"
          variant="outline"
          @click="resetToConfig"
        />
        <tp-button
          label="Done"
          @click="close"
        />
      </div>
    </div>
  </tp-dialog>
</template>

<script setup lang="ts">
import { ref, watch } from 'vue';
import TpButton from 'components/TpButton.vue';
import TpButtonGroup from 'components/TpButtonGroup.vue';
import TpDialog from 'components/TpDialog.vue';
import TpGlassInput from 'components/TpGlassInput.vue';
import TpIcon from 'components/TpIcon.vue';
import { api } from 'src/boot/axios';

interface ShareData {
  share_token: string;
  password: string;
  expires_at: string;
}

const props = defineProps<{
  modelValue?: boolean;
  jobId: string;
}>();

const emit = defineEmits<{
  (e: 'update:modelValue', value: boolean): void;
}>();

const loading = ref(false);
const error = ref<string | null>(null);
const shareData = ref<ShareData | null>(null);
const copiedLink = ref(false);
const copiedPassword = ref(false);

const shareUrl = ref('');

// Password configuration
const passwordMode = ref<'generate' | 'custom'>('generate');
const customPassword = ref('');

// Reset to initial state when dialog opens
watch(() => props.modelValue, (isOpen) => {
  if (isOpen) {
    // Reset all state when opening
    resetToConfig();
  }
});

function resetToConfig() {
  shareData.value = null;
  error.value = null;
  loading.value = false;
  copiedLink.value = false;
  copiedPassword.value = false;
  // Keep password mode and custom password so user can retry with same settings
}

async function generateShareLink() {
  loading.value = true;
  error.value = null;
  shareData.value = null;

  try {
    const payload = passwordMode.value === 'custom' && customPassword.value
      ? { custom_password: customPassword.value }
      : {};
    
    const response = await api.post(`/jobs/share/${props.jobId}`, payload);
    shareData.value = response.data;
    
    // Build the share URL
    const baseUrl = window.location.origin;
    shareUrl.value = `${baseUrl}/shared/${response.data.share_token}`;
  } catch (err: unknown) {
    const apiError = err as { response?: { status?: number; data?: { detail?: string } } };
    if (apiError.response?.status === 404) {
      error.value = 'Job not found. It may have been deleted.';
    } else {
      error.value = 'Failed to generate share link. Please try again.';
    }
    console.error('Error generating share link:', err);
  } finally {
    loading.value = false;
  }
}

function close() {
  emit('update:modelValue', false);
  // Reset state after a delay for smooth animation
  setTimeout(() => {
    shareData.value = null;
    error.value = null;
    copiedLink.value = false;
    copiedPassword.value = false;
    customPassword.value = '';
    passwordMode.value = 'generate';
  }, 300);
}

async function copyLink() {
  try {
    await navigator.clipboard.writeText(shareUrl.value);
    copiedLink.value = true;
    setTimeout(() => copiedLink.value = false, 2000);
  } catch (err) {
    console.error('Failed to copy link:', err);
  }
}

async function copyPassword() {
  if (!shareData.value) return;
  try {
    await navigator.clipboard.writeText(shareData.value.password);
    copiedPassword.value = true;
    setTimeout(() => copiedPassword.value = false, 2000);
  } catch (err) {
    console.error('Failed to copy password:', err);
  }
}

function formatExpiry(isoDate: string): string {
  const date = new Date(isoDate);
  return date.toLocaleDateString(undefined, {
    year: 'numeric',
    month: 'long',
    day: 'numeric',
  });
}
</script>

<style scoped lang="scss">
@use 'src/css/helpers/glass' as *;

.tp-description {
  font-size: 14px;
  color: var(--text-medium);
  line-height: 1.5;
  margin-bottom: 20px;
}

.tp-custom-password {
  display: flex;
  flex-direction: column;
  gap: 8px;
  margin-bottom: 20px;
}

.tp-generate-info {
  margin-bottom: 20px;
}

.tp-info-glass {
  @include glass-info-box;

  :deep(svg) {
    color: var(--text-brand-regular);
    flex-shrink: 0;
  }
}

.tp-share-success {
  display: flex;
  flex-direction: column;
  gap: 20px;
}

.tp-success-badge {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 10px;
  padding: 12px 20px;
  background: var(--surface-brand-extra-light);
  border-radius: 8px;
  color: var(--text-brand-regular);
  font-weight: 500;
  font-size: 14px;
}

.tp-badge-icon {
  display: flex;
  align-items: center;
  justify-content: center;
}

.tp-share-field {
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.tp-label {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 14px;
  font-weight: 500;
  color: var(--text);
}

.tp-copy-field {
  display: flex;
  gap: 8px;
  align-items: stretch;
}

.tp-copy-btn {
  @include glass;
  display: flex;
  align-items: center;
  gap: 6px;
  padding: 12px 16px;
  border-radius: 10px;
  color: var(--text);
  font-size: 13px;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.2s ease;
  white-space: nowrap;

  &:hover {
    background: rgba(255, 255, 255, 0.9);
    border-color: var(--stroke-light);
  }

  &.copied {
    background: var(--surface-brand-extra-light);
    border-color: var(--stroke-brand-regular);
    color: var(--text-brand-regular);
  }
}

.tp-hint {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
  color: var(--text-medium);
}

.tp-expiry-notice {
  @include glass-info-box;
}

.tp-actions {
  display: flex;
  justify-content: flex-end;
  gap: 12px;
  margin-top: 4px;
}
</style>
