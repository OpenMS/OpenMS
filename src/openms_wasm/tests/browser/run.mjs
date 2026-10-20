// Browser tests of the OpenMS WebAssembly module in Chromium (driven with Playwright).
//
//   node tests/browser/run.mjs <package directory> [--benchmark]
//
// The package directory holds openms.mjs (the package of one build: @openms/openms-wasm or
// @openms/openms-wasm-mt) or index.mjs (the assembled package @openms/openms-wasm with both
// builds: its loader has to start the threaded build on isolated pages, else the
// single-threaded one, with a warning when more than one thread was requested).
// Two local web servers serve it with the test pages and data: one with the headers that
// make a page cross-origin isolated (COOP/COEP; shared memory, hence threads, available),
// one without. Each scenario loads the module on the page or in a Web Worker and checks:
// that it starts (a threaded module with a fixed worker pool would hang on a page that is
// not isolated), the threads it reports, and the results of feature finding and a database
// search; some also search two files with searchFiles() file by file (provideFile, in a Web Worker with
// FileReaderSync), which has to give the results of the search with both files there, or 12 files at once.
// --benchmark adds a larger search at 1, 2 and 4 threads and prints the timings.
//
// Environment: PLAYWRIGHT_MODULE  the playwright or playwright-core package (path or name;
//                                 default: playwright)
//              CHROMIUM_PATH      the browser executable (default: the one of Playwright)

import fs from "node:fs";
import http from "node:http";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const sourceRoot = path.resolve(here, "../../../..");
const args = process.argv.slice(2);
const packageDir = path.resolve(args.find((a) => !a.startsWith("--")) ?? "");
const benchmark = args.includes("--benchmark");
const assembled = fs.existsSync(path.join(packageDir, "index.mjs"));
if (!assembled && !fs.existsSync(path.join(packageDir, "openms.mjs"))) {
  console.error("usage: node run.mjs <directory of openms.mjs or of the assembled package> [--benchmark]");
  process.exit(2);
}
const packageName = JSON.parse(fs.readFileSync(path.join(packageDir, "package.json"), "utf8")).name;
const threadedBuild = !assembled && packageName.endsWith("-mt");
const withThreads = assembled || threadedBuild;

//------------------------------------------------------------------------------
// Web servers

const roots = {
  "/pkg/": packageDir,
  "/test/": here,
  "/data/topp/": path.join(sourceRoot, "src/tests/topp"),
  "/data/examples/": path.join(sourceRoot, "share/OpenMS/examples/TOPPAS/data"),
};
const types = { ".html": "text/html", ".mjs": "text/javascript", ".js": "text/javascript", ".wasm": "application/wasm" };

function serve(isolated) {
  const server = http.createServer((request, response) => {
    const url = decodeURIComponent(new URL(request.url, "http://localhost").pathname);
    const prefix = Object.keys(roots).find((p) => url.startsWith(p));
    const file = prefix && path.resolve(roots[prefix], url.slice(prefix.length));
    if (!file || !file.startsWith(roots[prefix]) || !fs.existsSync(file) || !fs.statSync(file).isFile()) {
      response.writeHead(404).end();
      return;
    }
    const headers = { "Content-Type": types[path.extname(file)] ?? "application/octet-stream", "Cache-Control": "no-store" };
    if (isolated) {
      headers["Cross-Origin-Opener-Policy"] = "same-origin";
      headers["Cross-Origin-Embedder-Policy"] = "require-corp";
    }
    response.writeHead(200, headers);
    fs.createReadStream(file).pipe(response);
  });
  return new Promise((resolve) => server.listen(0, "127.0.0.1", () => resolve(server)));
}

//------------------------------------------------------------------------------
// Scenarios and what they have to give

// search.test.mjs: TOPP_SimpleSearchEngine_1
const expectedPeptides = ["DFASSGGYVLHLHR", "IALSRPNVEVVALNDPFITNDYAAYM(Oxidation)FK", "RPGADSDIGGFGGLFDLAQAGFR"];

const scenarios = [
  { isolated: false, host: "worker", numThreads: 4, files: true },
  { isolated: false, host: "window" },
  { isolated: true, host: "worker", files: true },
  { isolated: true, host: "worker", numThreads: 1 },
  { isolated: true, host: "window", numThreads: 2, files: true },
];
if (withThreads) scenarios.push({ isolated: true, host: "worker", numThreads: 4, manyFiles: 12 });
if (benchmark) {
  for (const n of withThreads ? [1, 2, 4] : [1]) {
    scenarios.push({ isolated: true, host: "worker", numThreads: n, benchmark: true });
  }
}

/// The threads a scenario has to report.
function expectedThreads(s, cores) {
  if (!withThreads || (assembled && !s.isolated)) return { threaded: false, available: false, numThreads: 1 };
  if (!s.isolated) return { threaded: true, available: false, numThreads: 1 };
  return { threaded: true, available: true, numThreads: s.numThreads ?? Math.min(4, Math.ceil(cores / 2)) };
}

/// The warning the loader has to give: more than one thread requested where there are none.
function expectedWarnings(s) {
  return assembled && !s.isolated && s.numThreads > 1 ? 1 : 0;
}

//------------------------------------------------------------------------------

const require = createRequire(import.meta.url);
const { chromium } = require(process.env.PLAYWRIGHT_MODULE ?? "playwright");
const browser = await chromium.launch({ executablePath: process.env.CHROMIUM_PATH || undefined });
const servers = { true: await serve(true), false: await serve(false) };
console.log(`${packageName}${assembled ? " (assembled, loader)" : ""} in ${browser.browserType().name()} ${browser.version()}`);

let failures = 0;
let reference;
for (const s of scenarios) {
  const name = `${s.isolated ? "isolated" : "not isolated"}, ${s.host}, numThreads ${s.numThreads ?? "default"}${s.files ? ", files" : ""}` +
    `${s.manyFiles ? `, ${s.manyFiles} files` : ""}${s.benchmark ? ", benchmark" : ""}`;
  const problems = [];
  const context = await browser.newContext();
  const page = await context.newPage();
  page.on("pageerror", (e) => problems.push(`page error: ${e.message}`));
  const query = new URLSearchParams({ module: assembled ? "/pkg/index.mjs" : "/pkg/openms.mjs", host: s.host });
  if (s.numThreads) query.set("numThreads", s.numThreads);
  if (s.benchmark) query.set("benchmark", "1");
  if (s.files) query.set("fileByFile", "1");
  if (s.manyFiles) query.set("manyFiles", s.manyFiles);
  const port = servers[s.isolated].address().port;
  let outcome;
  try {
    await page.goto(`http://127.0.0.1:${port}/test/page.html?${query}`);
    await page.waitForFunction(() => window.testResult, null, { timeout: s.benchmark ? 600000 : 180000, polling: 250 });
    outcome = await page.evaluate(() => window.testResult);
  } catch (e) {
    problems.push(`no result: ${e.message.split("\n")[0]} (a hang?)`);
  }
  await context.close();

  if (outcome?.error) problems.push(outcome.error);
  const r = outcome?.result;
  if (r) {
    const expected = expectedThreads(s, outcome.page.hardwareConcurrency);
    for (const key of Object.keys(expected)) {
      if (r.info[key] !== expected[key]) problems.push(`getThreadingInfo().${key} is ${r.info[key]}, expected ${expected[key]}`);
    }
    if (outcome.page.crossOriginIsolated !== s.isolated) problems.push(`crossOriginIsolated is ${outcome.page.crossOriginIsolated}`);
    if (r.warnings.length !== expectedWarnings(s) || r.warnings.some((w) => !w.includes("not cross-origin isolated"))) {
      problems.push(`warnings: ${JSON.stringify(r.warnings)}`);
    }
    if (r.featureTable.mz.length !== 13) problems.push(`${r.featureTable.mz.length} features, expected 13`);
    reference ??= r.featureTable;
    if (JSON.stringify(r.featureTable) !== JSON.stringify(reference)) problems.push("features differ from the first scenario");
    const peptides = r.search.hits.map((h) => h.sequence);
    if (r.search.status !== "EXECUTION_OK" || JSON.stringify(peptides) !== JSON.stringify(expectedPeptides)) {
      problems.push(`search: ${r.search.status}, ${JSON.stringify(peptides)}`);
    }
    if (r.prose.status !== "EXECUTION_OK" || JSON.stringify(r.prose.hits) !== JSON.stringify(expectedPeptides)) {
      problems.push(`ProSE: ${r.prose.status}, ${JSON.stringify(r.prose.hits)}`);
    }
    if (!r.prose.progressLabels.includes("Building fragment index...")) {
      problems.push(`ProSE progress: ${JSON.stringify(r.prose.progressLabels)}`);
    }
    if (s.files) {
      const f = r.files;
      if (!f || f.psms === 0 || !f.equal) problems.push(`searchFiles file by file: ${JSON.stringify(f)} (results equal to those with all files there: ${f?.equal})`);
      if (f && JSON.stringify(f.calls) !== JSON.stringify(["provide 0", "release 0", "provide 1", "release 1"])) problems.push(`searchFiles file by file: calls ${JSON.stringify(f.calls)}`);
      if (f && f.left.length) problems.push(`searchFiles file by file: files left ${JSON.stringify(f.left)}`);
      if (f && f.fileReaderSync !== (s.host === "worker")) problems.push(`searchFiles file by file: FileReaderSync ${f.fileReaderSync} on the ${s.host}`);
    }
    if (s.manyFiles && (r.manyFiles?.files !== s.manyFiles || !(r.manyFiles.psms > 0))) {
      problems.push(`searchFiles of ${s.manyFiles} files: ${JSON.stringify(r.manyFiles)}`);
    }
    // in a Web Worker, OpenMS must not block the page (the window host blocks it by design)
    if (s.host === "worker" && outcome.maxStallMs > 500) problems.push(`page blocked for ${outcome.maxStallMs} ms`);
  }

  failures += problems.length ? 1 : 0;
  const t = r?.timings ?? {};
  const details = r ? `threads ${r.info.numThreads}; start ${t.start} ms, features ${t.featureFinding} ms, search ${t.search} ms, ProSE ${t.prose} ms` +
    (r.files ? `, searchFiles file by file ${t.files} ms (${r.files.psms} PSMs)` : "") +
    (r.manyFiles ? `, searchFiles of ${r.manyFiles.files} files ${t.manyFiles} ms (${r.manyFiles.psms} PSMs)` : "") +
    (r.benchmark ? `, benchmark search ${t.benchmarkSearch} ms (${r.benchmark.psms} PSMs)` : "") +
    `; page blocked at most ${outcome.maxStallMs} ms` : "";
  console.log(`${problems.length ? "FAIL" : "ok  "} ${name}: ${details}`);
  for (const p of problems) console.log(`       ${p}`);
}

// The demo application (examples/browser, not in the assembled package) with its service worker
// (examples/sw.js) on the server without isolation headers: the worker has to serve the page, make
// it cross-origin isolated (threads), and keep serving it on a second visit (the module from its store).
let total = scenarios.length;
if (fs.existsSync(path.join(packageDir, "examples/browser/index.html"))) {
  ++total;
  const port = servers[false].address().port;
  const context = await browser.newContext();
  const problems = [];
  let status = "";
  for (let visit = 1; visit <= 2; ++visit) {
    const page = await context.newPage();
    page.on("pageerror", (e) => problems.push(`visit ${visit}: page error: ${e.message}`));
    try {
      await page.goto(`http://127.0.0.1:${port}/pkg/examples/browser/index.html`); // (the server serves files only)
      // (the page reloads itself once on the first visit: poll instead of waiting on one document)
      const deadline = Date.now() + 180000;
      let state = null;
      while (Date.now() < deadline) {
        try {
          state = await page.evaluate(() => ({
            status: document.getElementById("status").textContent,
            isolated: crossOriginIsolated,
            controlled: !!navigator.serviceWorker.controller,
          }));
          if (/is ready|could not|failed/i.test(state.status)) break;
        } catch {
          // navigating
        }
        await page.waitForTimeout(250);
      }
      status = state?.status ?? "(no status)";
      if (!state?.controlled) problems.push(`visit ${visit}: the service worker does not serve the page`);
      if (!state?.isolated) problems.push(`visit ${visit}: the page is not cross-origin isolated`);
      if (!/is ready/.test(status)) problems.push(`visit ${visit}: ${status}`);
    } catch (e) {
      problems.push(`visit ${visit}: ${e.message.split("\n")[0]}`);
    }
    await page.close();
  }
  await context.close();
  failures += problems.length ? 1 : 0;
  console.log(`${problems.length ? "FAIL" : "ok  "} demo application with its service worker, not isolated server: ${status}`);
  for (const p of problems) console.log(`       ${p}`);
}

await browser.close();
for (const server of Object.values(servers)) server.close();
console.log(failures ? `${failures} of ${total} scenarios failed` : `all ${total} scenarios passed`);
process.exit(failures ? 1 : 0);
