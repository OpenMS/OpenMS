// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// Round-trip test for the HTML INIFileEditor (share/OpenMS/HTML/INIFileEditor.html).
//
// The HTML file contains a self-contained "PARAM CORE" script block (parser +
// serializer for Param XML) that runs both in browsers and under Node. This
// script extracts that block, then feeds every .ini file of the repository
// through parse() -> serialize() and checks
//
//   1. byte-identity for files that are already in the canonical format of the
//      current C++ writer (ParamXMLFile, schema 1.8.0) -- i.e. the HTML editor
//      saves exactly the bytes ParamXMLFile::store() would produce;
//   2. idempotence for ALL parseable files: a second parse+serialize pass over
//      the editor's own output must reproduce it byte-for-byte (older files
//      are canonicalized once, exactly like a load+store through the C++ code
//      would canonicalize them).
//
// Usage: node tools/test_INIFileEditor_html.mjs [--verbose]

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const repo = join(dirname(fileURLToPath(import.meta.url)), '..');
const verbose = process.argv.includes('--verbose');

// ---- extract and evaluate the PARAM CORE block from the HTML file ----------
const html = readFileSync(join(repo, 'share/OpenMS/HTML/INIFileEditor.html'), 'utf8');
const begin = html.indexOf('===== BEGIN OPENMS PARAM CORE =====');
const end = html.indexOf('===== END OPENMS PARAM CORE =====');
if (begin < 0 || end < 0) throw new Error('PARAM CORE markers not found in INIFileEditor.html');
const core = html.slice(html.indexOf('*/', begin) + 2, html.lastIndexOf('/*', end));
const ParamCore = new Function(`${core}; return ParamCore;`)();

// ---- collect .ini files ----------------------------------------------------
function* walk(dir) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    const st = statSync(p);
    if (st.isDirectory()) yield* walk(p);
    else if (name.endsWith('.ini')) yield p;
  }
}
const roots = ['src/tests', 'share/OpenMS'].map(d => join(repo, d));
const files = roots.flatMap(r => { try { return [...walk(r)]; } catch { return []; } });

// ---- run -------------------------------------------------------------------
let nFiles = 0, nParsed = 0, nByteIdentical = 0, nCanonicalized = 0;
let nNotParam = 0, nParseFail = 0, nIdempotenceFail = 0;
const parseFailures = [], idempotenceFailures = [], canonicalized = [];

function firstDiff(a, b) {
  const la = a.split('\n'), lb = b.split('\n');
  for (let i = 0; i < Math.max(la.length, lb.length); i++) {
    if (la[i] !== lb[i]) return { line: i + 1, a: la[i], b: lb[i] };
  }
  return null;
}

for (const f of files) {
  nFiles++;
  const rel = relative(repo, f);
  const bytes = new Uint8Array(readFileSync(f));
  const text = ParamCore.decodeBytes(bytes);
  let parsed;
  try {
    parsed = ParamCore.parseParam(text);
  } catch (err) {
    // files that are not Param XML at all (other formats also use ".ini")
    if (/root element is|expected <PARAMETERS>|empty document/.test(err.message)) { nNotParam++; continue; }
    nParseFail++;
    parseFailures.push({ rel, err: err.message });
    continue;
  }
  nParsed++;
  const out1 = ParamCore.serializeParam(parsed.root);

  // 1. byte identity against the original
  if (Buffer.compare(Buffer.from(ParamCore.encodeBytes(out1)), Buffer.from(bytes)) === 0) {
    nByteIdentical++;
  } else {
    nCanonicalized++;
    canonicalized.push({ rel, version: parsed.version, diff: firstDiff(text, out1) });
  }

  // 2. idempotence: our own output must round-trip byte-for-byte
  try {
    const out2 = ParamCore.serializeParam(ParamCore.parseParam(out1).root);
    if (out2 !== out1) {
      nIdempotenceFail++;
      idempotenceFailures.push({ rel, diff: firstDiff(out1, out2) });
    }
  } catch (err) {
    nIdempotenceFail++;
    idempotenceFailures.push({ rel, err: err.message });
  }
}

// ---- report ----------------------------------------------------------------
console.log(`scanned .ini files:              ${nFiles}`);
console.log(`  not Param XML (skipped):       ${nNotParam}`);
console.log(`  parsed:                        ${nParsed}`);
console.log(`  parse failures:                ${nParseFail}`);
console.log(`  byte-identical round trip:     ${nByteIdentical}`);
console.log(`  canonicalized (diff to input): ${nCanonicalized}`);
console.log(`  idempotence failures:          ${nIdempotenceFail}`);

if (parseFailures.length) {
  console.log('\n--- parse failures ---');
  for (const p of parseFailures.slice(0, verbose ? Infinity : 10)) console.log(`  ${p.rel}\n    ${p.err}`);
}
if (idempotenceFailures.length) {
  console.log('\n--- idempotence failures (BUGS) ---');
  for (const p of idempotenceFailures.slice(0, verbose ? Infinity : 10)) {
    console.log(`  ${p.rel}`);
    if (p.err) console.log(`    ${p.err}`);
    else console.log(`    line ${p.diff.line}\n      pass1: ${p.diff.a}\n      pass2: ${p.diff.b}`);
  }
}
if (canonicalized.length && verbose) {
  console.log('\n--- canonicalized files (first differing line each) ---');
  for (const c of canonicalized) {
    console.log(`  ${c.rel} (file version ${c.version})`);
    if (c.diff) console.log(`    line ${c.diff.line}\n      file:   ${JSON.stringify(c.diff.a)}\n      editor: ${JSON.stringify(c.diff.b)}`);
  }
}

process.exitCode = (nParseFail || nIdempotenceFail) ? 1 : 0;
