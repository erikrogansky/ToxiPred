/**
 * Pre-build ketcher-react + ketcher-standalone into a standalone IIFE bundle.
 * This avoids Rollup's CJS/ESM interop issues with these complex packages
 * (Web Workers, WASM, babel runtime helpers, etc.).
 *
 * The output is placed in public/ketcher/ and loaded via <script> in index.html.
 * The main app accesses ketcher via window.__ketcher globals.
 */
import { build } from 'esbuild';
import { mkdirSync, copyFileSync } from 'fs';
import { dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const outdir = `${__dirname}/../public/ketcher`;
mkdirSync(outdir, { recursive: true });

await build({
  entryPoints: [`${__dirname}/ketcher-entry.js`],
  bundle: true,
  format: 'iife',
  globalName: '__ketcher',
  outfile: `${outdir}/ketcher-bundle.js`,
  define: {
    'process.env.NODE_ENV': '"production"',
    'process.env': JSON.stringify({ NODE_ENV: 'production' }),
    'process.version': '"v20.0.0"',
    'process.platform': '"browser"',
    global: 'globalThis',
  },
  banner: {
    js: 'if(typeof process==="undefined"){var process={env:{NODE_ENV:"production"},version:"v20.0.0",platform:"browser",nextTick:function(f){Promise.resolve().then(f)}}}',
  },
  loader: {
    '.wasm': 'binary',
  },
  target: 'es2022',
  minify: true,
  sourcemap: false,
  logLevel: 'info',
});

// Overwrite with the full ketcher-react CSS (esbuild only extracts a small subset)
const ketcherCss = `${__dirname}/../node_modules/ketcher-react/dist/index.css`;
copyFileSync(ketcherCss, `${outdir}/ketcher-bundle.css`);

console.log('Ketcher bundle built successfully.');
