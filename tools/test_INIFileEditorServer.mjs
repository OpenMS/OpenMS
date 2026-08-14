// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// API + security test for share/OpenMS/SCRIPTS/INIFileEditorServer.py, the
// local helper that gives the HTML INIFileEditor in-place saving in browsers
// without the File System Access API (Safari, Firefox).
//
// Needs only Node >= 18 and python3. Usage: node tools/test_INIFileEditorServer.mjs

import { spawn } from 'node:child_process';
import http from 'node:http';
import { copyFileSync, readFileSync, mkdtempSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { tmpdir } from 'node:os';

const repo = join(dirname(fileURLToPath(import.meta.url)), '..');
let failed = 0;
const check = (name, ok) => { console.log((ok ? 'PASS' : 'FAIL') + '  ' + name); if (!ok) failed = 1; };

// raw request helper (unlike fetch, node:http sends any Host header we set)
function req(port, { method = 'GET', path = '/', headers = {}, body = null } = {}) {
  return new Promise((resolve, reject) => {
    if (body !== null) {
      body = Buffer.isBuffer(body) ? body : Buffer.from(body);
      headers = { 'Content-Length': body.length, ...headers };
    }
    const r = http.request({ host: '127.0.0.1', port, method, path, headers }, (res) => {
      const chunks = [];
      res.on('data', (c) => chunks.push(c));
      res.on('end', () => resolve({ status: res.statusCode, headers: res.headers, body: Buffer.concat(chunks) }));
    });
    r.on('error', reject);
    if (body !== null) r.write(body);
    r.end();
  });
}

// start the helper on a temp copy of a canonical INI
const dir = mkdtempSync(join(tmpdir(), 'ife-server-'));
const ini = join(dir, 'test.ini');
copyFileSync(join(repo, 'src/tests/topp/WRITE_INI_OUT.ini'), ini);
const server = spawn('python3', [join(repo, 'share/OpenMS/SCRIPTS/INIFileEditorServer.py'), '--no-browser', ini]);
const url = await new Promise((resolve, reject) => {
  let buf = '';
  server.stdout.on('data', (d) => { buf += d; const m = /Editor URL: (\S+)/.exec(buf); if (m) resolve(m[1]); });
  server.stderr.on('data', (d) => process.stderr.write(d));
  server.on('exit', (c) => reject(new Error('helper exited with code ' + c)));
  setTimeout(() => reject(new Error('no URL printed by helper')), 8000);
});
const port = Number(new URL(url).port);
const token = /#t=([A-Za-z0-9_-]+)/.exec(url)[1];
const auth = { 'X-Auth-Token': token };
const host = `127.0.0.1:${port}`;

try {
  // page + 404s
  check('serves editor page at /', (await req(port, { headers: { Host: host } })).status === 200);
  check('unknown paths are 404', (await req(port, { path: '/etc/passwd', headers: { Host: host } })).status === 404);
  check('POST has no endpoints', (await req(port, { method: 'POST', path: '/api/file/0', headers: { Host: host, ...auth } })).status === 404);

  // auth
  check('listing without token rejected', (await req(port, { path: '/api/files', headers: { Host: host } })).status === 403);
  check('listing with wrong token rejected', (await req(port, { path: '/api/files', headers: { Host: host, 'X-Auth-Token': 'nope' } })).status === 403);
  check('PUT without token rejected', (await req(port, { method: 'PUT', path: '/api/file/0', headers: { Host: host }, body: 'x' })).status === 403);

  // DNS rebinding / CSRF hardening
  check('spoofed Host rejected', (await req(port, { path: '/api/files', headers: { Host: 'evil.example', ...auth } })).status === 403);
  check('foreign Origin rejected', (await req(port, { path: '/api/files', headers: { Host: host, Origin: 'http://evil.example', ...auth } })).status === 403);
  check('own Origin accepted', (await req(port, { path: '/api/files', headers: { Host: host, Origin: `http://${host}`, ...auth } })).status === 200);

  // listing + read
  const listing = JSON.parse((await req(port, { path: '/api/files', headers: { Host: host, ...auth } })).body);
  check('listing contains the file', listing.files.length === 1 && listing.files[0].name === 'test.ini');
  const got = await req(port, { path: '/api/file/0', headers: { Host: host, ...auth } });
  check('file readable via API', got.status === 200 && got.body.equals(readFileSync(ini)));
  const mtime = got.headers['x-file-mtime'];
  check('mtime reported', typeof mtime === 'string' && mtime.length > 0);
  check('out-of-range file id is 404', (await req(port, { path: '/api/file/1', headers: { Host: host, ...auth } })).status === 404);

  // write: stale mtime -> 409, correct mtime -> 200, content lands on disk
  const newContent = got.body.toString('latin1').replace('value="4"', 'value="7"');
  const stale = await req(port, { method: 'PUT', path: '/api/file/0', headers: { Host: host, ...auth, 'X-Expected-Mtime': '1' }, body: newContent });
  check('stale mtime rejected with 409', stale.status === 409);
  const put = await req(port, { method: 'PUT', path: '/api/file/0', headers: { Host: host, ...auth, 'X-Expected-Mtime': mtime }, body: Buffer.from(newContent, 'latin1') });
  check('PUT with matching mtime succeeds', put.status === 200);
  check('file updated in place', readFileSync(ini, 'latin1') === newContent);
  check('new mtime returned', put.headers['x-file-mtime'] !== mtime);
  const forced = await req(port, { method: 'PUT', path: '/api/file/0', headers: { Host: host, ...auth }, body: Buffer.from(newContent, 'latin1') });
  check('PUT without expected mtime overwrites (forced save)', forced.status === 200);
} finally {
  server.kill('SIGINT');
}
process.exitCode = failed;
