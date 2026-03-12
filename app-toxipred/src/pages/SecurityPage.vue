<template>
  <tp-page class="tp-legal-page">
    <h1>Security</h1>
    <p class="paragraph tp-legal-page__updated">Last updated: {{ lastUpdated }}</p>

    <section>
      <h2>1. Overview</h2>
      <p class="paragraph">
        ToxiPred is committed to maintaining a secure platform. This page describes the security measures we
        implement and provides guidance for responsible disclosure if you discover a vulnerability.
      </p>
    </section>

    <section>
      <h2>2. Architecture &amp; Infrastructure</h2>

      <h3>2.1 Network Security</h3>
      <ul>
        <li><strong>HTTPS/TLS encryption</strong> — all communication between your browser and our servers is
          encrypted in transit using TLS. Plain HTTP connections are redirected to HTTPS.</li>
        <li><strong>Reverse-proxy isolation</strong> — a nginx reverse proxy sits between the public internet
          and our application services, providing an additional security layer, request filtering, and rate
          limiting.</li>
        <li><strong>No direct database exposure</strong> — the database is only accessible from internal
          application services and is not exposed to the public internet.</li>
      </ul>

      <h3>2.2 Application Security</h3>
      <ul>
        <li><strong>No authentication surface</strong> — ToxiPred does not have user accounts, login forms,
          or session tokens, which eliminates entire categories of attacks (credential stuffing, session
          hijacking, password breaches).</li>
        <li><strong>Input validation</strong> — all chemical identifiers submitted by users are validated and
          sanitized before processing. SMILES strings are parsed and canonicalized using RDKit, a
          well-established cheminformatics library.</li>
        <li><strong>No user-generated content rendering</strong> — the platform does not render arbitrary
          user-supplied HTML, Markdown, or other rich content, minimizing cross-site scripting (XSS) risk.</li>
        <li><strong>CORS policy</strong> — Cross-Origin Resource Sharing headers are configured to control
          which origins can interact with our API.</li>
      </ul>

      <h3>2.3 Data Security</h3>
      <ul>
        <li><strong>Minimal data collection</strong> — we collect the absolute minimum data necessary to
          provide predictions. No personal identifiers are stored. See our
          <router-link to="/privacy-policy">Privacy Policy</router-link> for details.</li>
        <li><strong>Automatic data deletion</strong> — prediction results are automatically purged from the
          server after 14 days via a scheduled background task. Share links expire after 30 days.</li>
        <li><strong>Password hashing</strong> — when you password-protect a shared prediction link, the
          password is cryptographically hashed before storage and never stored in plain text.</li>
        <li><strong>Cryptographic share tokens</strong> — share links use cryptographically secure random
          tokens that are infeasible to guess or enumerate.</li>
      </ul>
    </section>

    <section>
      <h2>3. Third-Party Dependencies</h2>
      <p class="paragraph">
        We regularly review and update our software dependencies to address known vulnerabilities. Key
        components include:
      </p>
      <ul>
        <li><strong>Backend:</strong> Python, FastAPI, SQLAlchemy, Celery, RDKit — maintained with
          security patches from upstream projects.</li>
        <li><strong>Frontend:</strong> Vue.js, Quasar Framework, TypeScript — bundled and served from
          our own infrastructure with no external CDN dependencies.</li>
        <li><strong>Infrastructure:</strong> Docker containers with minimal base images, nginx for reverse
          proxying, Redis for task queueing.</li>
      </ul>
      <p class="paragraph">
        All third-party assets (Ketcher molecular editor, XSMILES renderer, fonts) are bundled locally
        and served from our servers. No external scripts are loaded at runtime.
      </p>
    </section>

    <section>
      <h2>4. What We Do <em>Not</em> Do</h2>
      <p class="paragraph">
        To reduce the attack surface and protect your privacy, ToxiPred intentionally avoids:
      </p>
      <ul>
        <li>Storing any personal identifiers (names, emails, IP addresses) in the application database.</li>
        <li>Using cookies, tracking pixels, analytics scripts, or browser fingerprinting.</li>
        <li>Loading scripts, stylesheets, or fonts from third-party CDNs.</li>
        <li>Implementing user accounts or authentication flows (no passwords to breach).</li>
        <li>Executing user-supplied code or rendering user-supplied HTML.</li>
      </ul>
    </section>

    <section>
      <h2>5. Responsible Disclosure</h2>
      <div class="tp-legal-page__highlight">
        <p class="paragraph">
          If you discover a security vulnerability in ToxiPred, we appreciate your help in disclosing it
          responsibly. Please <strong>do not</strong> open a public issue or disclose the vulnerability
          publicly before we have had a chance to address it.
        </p>
      </div>

      <h3>5.1 How to Report</h3>
      <p class="paragraph">Send a detailed report to:</p>
      <ul>
        <li><tp-link href="mailto:roganskyerik@gmail.com" label="roganskyerik@gmail.com" underline="hover" />
          — with the subject line: <strong>Security Vulnerability Report</strong></li>
      </ul>

      <h3>5.2 What to Include</h3>
      <ul>
        <li>A description of the vulnerability and its potential impact.</li>
        <li>Steps to reproduce the issue.</li>
        <li>Affected component(s) — frontend, backend API, infrastructure, etc.</li>
        <li>Any proof-of-concept code or screenshots, if applicable.</li>
        <li>Your preferred contact method for follow-up (optional).</li>
      </ul>

      <h3>5.3 Our Commitment</h3>
      <ul>
        <li>We will acknowledge receipt of your report within 48 hours.</li>
        <li>We will investigate and provide a timeline for a fix.</li>
        <li>We will notify you when the vulnerability has been resolved.</li>
        <li>We will not take legal action against security researchers who act in good faith and follow
          this responsible disclosure process.</li>
      </ul>
    </section>

    <section>
      <h2>6. Scope of Testing</h2>
      <p class="paragraph">
        If you are testing ToxiPred for security vulnerabilities, please stay within reasonable bounds:
      </p>
      <ul>
        <li><strong>Do</strong> test the publicly accessible web application and API.</li>
        <li><strong>Do not</strong> perform denial-of-service (DoS/DDoS) attacks.</li>
        <li><strong>Do not</strong> attempt to access, modify, or delete other users' data (shared predictions).</li>
        <li><strong>Do not</strong> use automated vulnerability scanners that generate excessive traffic.</li>
        <li><strong>Do not</strong> attempt to compromise the underlying server infrastructure, operating system, or network.</li>
      </ul>
    </section>

    <section>
      <h2>7. Limitations</h2>
      <p class="paragraph">
        ToxiPred is an academic research project, not a commercial service with a dedicated security team. While
        we implement industry-standard security practices and respond promptly to reported issues, we make no
        guarantees regarding the absolute security of the platform. Users should consider this when deciding
        what data to submit.
      </p>
      <p class="paragraph">
        For sensitive or proprietary chemical structures, we recommend using the platform's SMILES-based input
        rather than compound names, as SMILES strings do not inherently reveal the identity of proprietary
        compounds to anyone without specialized knowledge.
      </p>
    </section>

    <section>
      <h2>8. Contact</h2>
      <p class="paragraph">For security-related inquiries:</p>
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
    border-left: 4px solid var(--accent-blue-regular);
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
}
</style>
