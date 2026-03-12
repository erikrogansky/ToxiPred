<template>
  <tp-page class="tp-legal-page">
    <h1>Privacy Policy</h1>
    <p class="paragraph tp-legal-page__updated">Last updated: {{ lastUpdated }}</p>

    <section>
      <h2>1. Introduction</h2>
      <p class="paragraph">
        ToxiPred ("<strong>we</strong>," "<strong>us</strong>," or "<strong>our</strong>") is a free, open-access
        research platform for predicting dermatological toxicity of chemical compounds using QSAR (Quantitative
        Structure–Activity Relationship) models. It is developed and maintained at the Slovak University of
        Technology in Bratislava (STU).
      </p>
      <p class="paragraph">
        This Privacy Policy explains what information we collect, how we use it, and your rights regarding that
        information. We are committed to protecting your privacy and process only the minimum data necessary
        to provide our service.
      </p>
    </section>

    <section>
      <h2>2. Information We Collect</h2>

      <h3>2.1 Chemical Data You Submit</h3>
      <p class="paragraph">
        When you request a toxicity prediction, you provide one or more chemical identifiers (SMILES strings,
        CAS numbers, or trivial names). We process this input to generate prediction results, which include:
      </p>
      <ul>
        <li>Canonical SMILES representation</li>
        <li>Molecular formula and compound names</li>
        <li>Prediction outcomes (toxic / non-toxic classification)</li>
        <li>Confidence scores and applicability-domain assessment</li>
        <li>Feature importance and atom-level contribution scores</li>
      </ul>

      <h3>2.2 Technical Data</h3>
      <p class="paragraph">
        Our web server (nginx) may receive standard HTTP metadata such as your IP address, browser user-agent
        string, and request timestamps as part of normal network communication. This data is handled in transit
        by the reverse-proxy layer but is <strong>not stored, logged, or associated with your predictions</strong>
        in the application database.
      </p>

      <h3>2.3 Data We Do <em>Not</em> Collect</h3>
      <ul>
        <li>We do <strong>not</strong> require registration, accounts, or authentication of any kind.</li>
        <li>We do <strong>not</strong> collect email addresses, names, phone numbers, or any other personal identifiers.</li>
        <li>We do <strong>not</strong> use cookies.</li>
        <li>We do <strong>not</strong> use analytics, tracking pixels, fingerprinting, or any third-party tracking services.</li>
        <li>We do <strong>not</strong> share any data with third parties.</li>
      </ul>
    </section>

    <section>
      <h2>3. How We Use Your Data</h2>
      <p class="paragraph">Chemical data you submit is used solely to:</p>
      <ul>
        <li>Run QSAR toxicity prediction models and return results to you.</li>
        <li>Cache results server-side so that repeated queries for the same compound return faster.</li>
        <li>Enable the optional sharing feature if you choose to generate a share link.</li>
      </ul>
      <p class="paragraph">
        We do not perform data mining, profiling, or automated decision-making on your submissions. We do not
        use your chemical queries to train or improve our models.
      </p>
    </section>

    <section>
      <h2>4. Data Retention</h2>
      <table>
        <thead>
          <tr>
            <th>Data Type</th>
            <th>Storage Location</th>
            <th>Retention Period</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>Prediction results</td>
            <td>Server database</td>
            <td>14 days (auto-deleted)</td>
          </tr>
          <tr>
            <td>Share links</td>
            <td>Server database</td>
            <td>30 days (auto-expired)</td>
          </tr>
          <tr>
            <td>Job history</td>
            <td>Your browser (localStorage)</td>
            <td>14 days (auto-pruned on next visit)</td>
          </tr>
          <tr>
            <td>Accessibility preferences</td>
            <td>Your browser (localStorage)</td>
            <td>Until you clear browser data</td>
          </tr>
        </tbody>
      </table>
      <p class="paragraph">
        Server-side data is automatically and permanently deleted after the retention period by a scheduled
        background task. There is no manual override or recovery once data is deleted.
      </p>
    </section>

    <section>
      <h2>5. Local Storage</h2>
      <p class="paragraph">
        We use your browser's localStorage (not cookies) to store the following on your device:
      </p>
      <ul>
        <li><strong>Job history</strong> — a local record of your prediction results so you can revisit them in your Workspace without re-submitting.</li>
        <li><strong>Accessibility settings</strong> — your dark-mode, high-contrast, and font-family preferences.</li>
      </ul>
      <p class="paragraph">
        This data never leaves your device and is not transmitted to our servers. You can clear it at any time
        through your browser settings.
      </p>
    </section>

    <section>
      <h2>6. Data Sharing Features</h2>
      <p class="paragraph">
        ToxiPred lets you share prediction results with others via a unique, time-limited link. When you create
        a share link:
      </p>
      <ul>
        <li>A cryptographically random token is generated to identify the shared job.</li>
        <li>You may optionally set a password. Passwords are hashed before storage and never stored in plain text.</li>
        <li>Share links expire automatically after 30 days.</li>
        <li>Only the specific prediction results you choose to share are accessible via the link.</li>
      </ul>
      <p class="paragraph">
        You are responsible for distributing share links only to intended recipients. Anyone with a valid link
        (and password, if set) can view the shared prediction data.
      </p>
    </section>

    <section>
      <h2>7. Data Security</h2>
      <p class="paragraph">
        We implement reasonable technical measures to protect data in transit and at rest, including:
      </p>
      <ul>
        <li>HTTPS/TLS encryption for all client–server communication.</li>
        <li>Reverse-proxy isolation between the public internet and application services.</li>
        <li>Automatic data deletion after retention periods.</li>
        <li>No storage of personal identifiers that could be linked to chemical queries.</li>
      </ul>
      <p class="paragraph">
        However, no system is perfectly secure. Please see our
        <router-link to="/security">Security</router-link> page for more details.
      </p>
    </section>

    <section>
      <h2>8. Third-Party Services</h2>
      <p class="paragraph">
        ToxiPred does <strong>not</strong> integrate any third-party analytics, advertising, social-media, or
        tracking services. All computational resources (QSAR models, cheminformatics libraries, molecular
        editors) run on our own infrastructure or are bundled locally in the application.
      </p>
    </section>

    <section>
      <h2>9. Children's Privacy</h2>
      <p class="paragraph">
        ToxiPred is a scientific research tool and is not directed at children under 16. We do not knowingly
        collect any personal information from anyone, including children. Because no personal data is collected,
        no age-specific data handling applies.
      </p>
    </section>

    <section>
      <h2>10. Your Rights</h2>
      <p class="paragraph">
        Because ToxiPred does not collect personal data or use accounts, traditional data-subject rights
        (access, rectification, erasure, portability) under GDPR or similar regulations have limited
        applicability. Nonetheless:
      </p>
      <ul>
        <li>You can clear all locally stored data at any time by clearing your browser's localStorage.</li>
        <li>Server-side prediction results are automatically deleted after 14 days.</li>
        <li>If you have concerns about any data associated with your use, contact us and we will address them promptly.</li>
      </ul>
    </section>

    <section>
      <h2>11. Changes to This Policy</h2>
      <p class="paragraph">
        We may update this Privacy Policy from time to time. Changes will be reflected by updating the
        "Last updated" date at the top of this page. Continued use of ToxiPred after changes constitutes
        acceptance of the revised policy.
      </p>
    </section>

    <section>
      <h2>12. Contact</h2>
      <p class="paragraph">If you have questions about this Privacy Policy, contact us at:</p>
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
