<template>
  <q-dialog :model-value="modelValue ?? false" @update:model-value="emit('update:modelValue', $event)">
    <q-card class="tp-share-dialog">
      <q-card-section class="tp-dialog-header">
        <h3 class="tp-h3">Share Prediction</h3>
        <tp-icon-button
          icon-name="close-square"
          weight="regular"
          :size="24"
          @click="close"
        />
      </q-card-section>

      <q-card-section class="tp-dialog-content">
        <!-- Initial State: Configure Password -->
        <div v-if="!shareData && !loading && !error" class="tp-config-section">
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
            <div class="tp-password-input-wrapper">
              <input
                v-model="customPassword"
                :type="showCustomPassword ? 'text' : 'password'"
                class="tp-glass-input"
                placeholder="Enter a secure password..."
              />
              <button class="tp-visibility-btn" @click="showCustomPassword = !showCustomPassword">
                <tp-icon :icon-name="showCustomPassword ? 'eye-slash' : 'eye'" weight="regular" :size="18" />
              </button>
            </div>
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

        <!-- Loading State -->
        <div v-if="loading" class="tp-loading-state">
          <div class="tp-loading-spinner">
            <q-spinner-dots size="40px" color="primary" />
          </div>
          <p class="tp-loading-text">Creating share link...</p>
        </div>

        <!-- Error State -->
        <div v-else-if="error" class="tp-error-state">
          <div class="tp-error-icon">
            <tp-icon icon-name="warning-2" weight="regular" :size="48" />
          </div>
          <p class="tp-error-text">{{ error }}</p>
          <div class="tp-actions">
            <tp-button
              label="Try Again"
              variant="outline"
              @click="resetToConfig"
            />
          </div>
        </div>

        <!-- Success State -->
        <div v-else-if="shareData" class="tp-share-success">
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
              <input
                type="text"
                :value="shareUrl"
                readonly
                class="tp-glass-input"
                ref="linkInput"
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
              <input
                type="text"
                :value="shareData.password"
                readonly
                class="tp-glass-input tp-password-display"
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
      </q-card-section>
    </q-card>
  </q-dialog>
</template>

<script setup lang="ts">
import { ref, watch } from 'vue';
import TpButton from 'components/TpButton.vue';
import TpButtonGroup from 'components/TpButtonGroup.vue';
import TpIconButton from 'components/TpIconButton.vue';
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
const showCustomPassword = ref(false);

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
.tp-share-dialog {
  min-width: 460px;
  max-width: 520px;
  border-radius: 20px;
  background: var(--glass-background);
  backdrop-filter: var(--glass-blur);
  border: 1px solid var(--glass-border);
  box-shadow: var(--glass-shadow);
  
  @media (max-width: 500px) {
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

.tp-description {
  font-size: 14px;
  color: var(--text-medium);
  line-height: 1.5;
  margin-bottom: 20px;
}

// Custom Password Input
.tp-custom-password {
  display: flex;
  flex-direction: column;
  gap: 8px;
  margin-bottom: 20px;
}

.tp-password-input-wrapper {
  display: flex;
  align-items: center;
  gap: 8px;
  border: 1px solid var(--stroke-light);
  border-radius: 8px;
  padding-right: 8px;
  
  &:focus-within {
    border-color: var(--stroke-brand-regular);
  }
}

.tp-visibility-btn {
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 6px;
  border: none;
  background: transparent;
  color: var(--text-light);
  cursor: pointer;
  border-radius: 6px;
  transition: all 0.15s ease;
  
  &:hover {
    color: var(--text);
  }
}

// Generate Info Box
.tp-generate-info {
  margin-bottom: 20px;
}

.tp-info-glass {
  display: flex;
  align-items: center;
  gap: 12px;
  padding: 14px 16px;
  background: var(--glass-background);
  backdrop-filter: blur(8px);
  border: 1px solid var(--glass-border);
  border-radius: 12px;
  font-size: 13px;
  color: var(--text-medium);
  
  :deep(svg) {
    color: var(--text-brand-regular);
    flex-shrink: 0;
  }
}

// Loading State
.tp-loading-state {
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 32px 0;
}

.tp-loading-spinner {
  margin-bottom: 16px;
}

.tp-loading-text {
  font-size: 14px;
  color: var(--text-medium);
}

// Error State
.tp-error-state {
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 24px 0;
}

.tp-error-icon {
  margin-bottom: 16px;
  color: var(--text-negative);
}

.tp-error-text {
  font-size: 14px;
  color: var(--text-negative);
  text-align: center;
  margin-bottom: 20px;
}

// Success State
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

// Share Fields
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

.tp-glass-input {
  flex: 1;
  padding: 12px 16px;
  border: 1px solid var(--glass-border);
  border-radius: 10px;
  font-size: 14px;
  background: var(--glass-background);
  backdrop-filter: blur(8px);
  color: var(--text);
  transition: all 0.2s ease;
  
  &::placeholder {
    color: var(--text-medium);
  }
  
  &:focus {
    outline: none;
    border-color: var(--stroke-brand-regular);
    box-shadow: 0 0 0 3px rgba(var(--color-brand-rgb), 0.1);
  }
  
  &:read-only {
    cursor: default;
    background: var(--surface-gray-light);
  }
}

.tp-password-display {
  font-family: monospace;
  letter-spacing: 0.5px;
}

.tp-copy-btn {
  display: flex;
  align-items: center;
  gap: 6px;
  padding: 12px 16px;
  border: 1px solid var(--glass-border);
  border-radius: 10px;
  background: var(--glass-background);
  backdrop-filter: blur(8px);
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
  display: flex;
  align-items: center;
  gap: 10px;
  padding: 14px 16px;
  background: var(--glass-background);
  backdrop-filter: blur(8px);
  border: 1px solid var(--glass-border);
  border-radius: 12px;
  font-size: 13px;
  color: var(--text-medium);
}

.tp-actions {
  display: flex;
  justify-content: flex-end;
  gap: 12px;
  margin-top: 4px;
}
</style>
