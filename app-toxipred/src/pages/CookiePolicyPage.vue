<template>
  <tp-page class="tp-legal-page">
    <h1>Cookie Policy</h1>
    <p class="paragraph tp-legal-page__updated">Last updated: {{ lastUpdated }}</p>

    <section>
      <h2>1. The Short Version</h2>
      <div class="tp-legal-page__highlight">
        <p class="paragraph bold">ToxiPred does not use cookies.</p>
        <p class="paragraph">
          We do not set any cookies — no first-party cookies, no third-party cookies, no tracking
          cookies, no advertising cookies, and no analytics cookies. You will never see a cookie consent
          banner on ToxiPred because there is nothing to consent to.
        </p>
      </div>
    </section>

    <section>
      <h2>2. What We Use Instead</h2>
      <p class="paragraph">
        To provide a functional user experience, we use your browser's
        <strong>localStorage</strong> — a client-side storage mechanism that keeps data exclusively on your
        device. Unlike cookies, localStorage data is never sent to our servers with HTTP requests.
      </p>

      <h3>2.1 What Is Stored Locally</h3>
      <table>
        <thead>
          <tr>
            <th>Key</th>
            <th>Purpose</th>
            <th>Contains</th>
            <th>Lifetime</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><code>toxipred.jobs.v1</code></td>
            <td>Prediction history</td>
            <td>Job IDs, SMILES, compound names, formulas, model names, prediction results, timestamps</td>
            <td>14 days (auto-pruned)</td>
          </tr>
          <tr>
            <td><code>darkMode</code></td>
            <td>Theme preference</td>
            <td>Boolean (on/off)</td>
            <td>Until cleared</td>
          </tr>
          <tr>
            <td><code>highContrast</code></td>
            <td>High-contrast mode</td>
            <td>Boolean (on/off)</td>
            <td>Until cleared</td>
          </tr>
          <tr>
            <td><code>fontFamily</code></td>
            <td>Accessibility font</td>
            <td>Font name (e.g., "Inter," "OpenDyslexic")</td>
            <td>Until cleared</td>
          </tr>
        </tbody>
      </table>
    </section>

    <section>
      <h2>3. How localStorage Differs from Cookies</h2>
      <table>
        <thead>
          <tr>
            <th>Feature</th>
            <th>Cookies</th>
            <th>localStorage</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>Sent with every HTTP request</td>
            <td>Yes</td>
            <td>No</td>
          </tr>
          <tr>
            <td>Accessible by our server</td>
            <td>Yes</td>
            <td>No</td>
          </tr>
          <tr>
            <td>Can be used for cross-site tracking</td>
            <td>Yes (third-party)</td>
            <td>No</td>
          </tr>
          <tr>
            <td>Requires consent banner (under GDPR)</td>
            <td>Usually yes</td>
            <td>Not for essential functionality</td>
          </tr>
          <tr>
            <td>Data stays on your device</td>
            <td>Not always</td>
            <td>Always</td>
          </tr>
        </tbody>
      </table>
      <p class="paragraph">
        Because localStorage data never leaves your browser, it poses no privacy risk to you and cannot be
        used for tracking or profiling.
      </p>
    </section>

    <section>
      <h2>4. Third-Party Cookies</h2>
      <p class="paragraph">
        ToxiPred does not embed any third-party services that could set cookies. We do not use:
      </p>
      <ul>
        <li>Analytics platforms (Google Analytics, Matomo, etc.)</li>
        <li>Advertising networks</li>
        <li>Social media widgets or share buttons</li>
        <li>External CDNs that require cookie-based authentication</li>
        <li>A/B testing or personalization services</li>
      </ul>
      <p class="paragraph">
        All assets (JavaScript libraries, fonts, icons, the Ketcher molecular editor) are bundled and served
        directly from our own infrastructure.
      </p>
    </section>

    <section>
      <h2>5. How to Clear Local Data</h2>
      <p class="paragraph">
        You can remove all data ToxiPred stores on your device at any time:
      </p>
      <ul>
        <li><strong>Chrome:</strong> Settings → Privacy and Security → Clear Browsing Data → check "Cookies and other site data" (this also covers localStorage).</li>
        <li><strong>Firefox:</strong> Settings → Privacy &amp; Security → Cookies and Site Data → Clear Data.</li>
        <li><strong>Safari:</strong> Preferences → Privacy → Manage Website Data → search for ToxiPred → Remove.</li>
        <li><strong>Edge:</strong> Settings → Privacy, search, and services → Clear browsing data.</li>
      </ul>
      <p class="paragraph">
        Clearing this data will reset your Workspace (prediction history) and accessibility preferences to
        defaults. It will not affect data already stored on our server (which auto-deletes after 14 days).
      </p>
    </section>

    <section>
      <h2>6. GDPR &amp; ePrivacy Compliance</h2>
      <p class="paragraph">
        Under the EU ePrivacy Directive and GDPR, consent is required for non-essential cookies and similar
        tracking technologies. Because ToxiPred:
      </p>
      <ul>
        <li>Uses no cookies whatsoever,</li>
        <li>Uses localStorage only for essential functionality (persisting your own prediction history and
          accessibility preferences),</li>
        <li>Does not track, profile, or fingerprint users,</li>
        <li>Does not share data with third parties,</li>
      </ul>
      <p class="paragraph">
        we believe no cookie consent mechanism is required. The localStorage usage falls under the "strictly
        necessary" exemption as it is required for the service to function as requested by the user.
      </p>
    </section>

    <section>
      <h2>7. Changes to This Policy</h2>
      <p class="paragraph">
        If we ever introduce cookies or change our use of browser storage, we will update this page and, if
        required by law, implement an appropriate consent mechanism before doing so.
      </p>
    </section>

    <section>
      <h2>8. Contact</h2>
      <p class="paragraph">If you have questions about this Cookie Policy, contact us at:</p>
      <ul>
        <li><tp-link href="mailto:roganskyerik@gmail.com" label="roganskyerik@gmail.com" underline="hover" /></li>
        <li><tp-link href="mailto:marta.prnova@stuba.sk" label="marta.prnova@stuba.sk" underline="hover" /></li>
      </ul>
    </section>
  </tp-page>
</template>

<script setup lang="ts">
import TpPage from 'components/TpPage.vue';
import TpLink from 'components/TpLink.vue';

const lastUpdated = 'March 12, 2026';
</script>

<style scoped lang="scss">
@use 'src/css/helpers/mixins' as *;
@use 'src/css/helpers/glass' as glass;

.tp-legal-page {
  box-sizing: content-box;
  max-width: 750px;
  gap: 0;

  &__updated {
    color: var(--text-medium);
    margin-bottom: 32px;
  }

  &__highlight {
    @include glass.glass-info-box();
    border-left: 4px solid var(--surface-brand-regular);
    padding: 20px 24px;
    margin: 8px 0 16px;
  }

  section {
    margin-bottom: 32px;
  }

  h1 {
    margin-bottom: 8px;
  }

  h2 {
    margin-bottom: 12px;
    padding-top: 8px;
  }

  h3 {
    margin-top: 16px;
    margin-bottom: 8px;
  }

  p {
    margin-bottom: 12px;
    color: var(--text);
    line-height: 1.7;
  }

  ul {
    padding-left: 24px;
    margin-bottom: 12px;

    li {
      color: var(--text);
      line-height: 1.8;
      font-size: 16px;
    }
  }

  a {
    color: var(--text-brand-regular);
    text-decoration: none;

    &:hover {
      text-decoration: underline;
    }
  }

  code {
    background: var(--glass-background);
    padding: 2px 6px;
    border-radius: 4px;
    font-size: 13px;
    font-family: monospace;
  }

  table {
    width: 100%;
    border-collapse: collapse;
    margin: 16px 0;
    @include glass.glass($radius: 12px);
    overflow: hidden;

    @include down(sm) {
      display: block;
      overflow-x: auto;
      -webkit-overflow-scrolling: touch;
    }

    th, td {
      padding: 12px 16px;
      text-align: left;
      font-size: 14px;
      color: var(--text);
      border-bottom: 1px solid var(--glass-border);
      white-space: nowrap;

      @include down(sm) {
        padding: 10px 12px;
        font-size: 13px;
      }
    }

    th {
      font-weight: 700;
      background: var(--glass-background);
    }

    tr:last-child td {
      border-bottom: none;
    }
  }
}
</style>
