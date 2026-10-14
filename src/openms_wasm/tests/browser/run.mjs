// Browser tests of the OpenMS WebAssembly module in Chromium (driven with Playwright).
//
//   node tests/browser/run.mjs <package directory> [--benchmark]
//
// The package directory holds openms.mjs (@openms/openms-wasm or @openms/openms-wasm-mt).
// Two local web servers serve it with the test pages and data: one with the headers that
// make a page cross-origin isolated (COOP/COEP; shared memory, hence threads, available),
// one without. Each scenario loads the module on the page or in a Web Worker and checks:
// that it starts (a threaded module with a fixed worker pool would hang on a page that is
// not isolated), the threads it reports, and the results of feature finding and a database
// search. --benchmark adds a larger search at 1, 2 and 4 threads and prints the timings.
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
if (!fs.existsSync(path.join(packageDir, "openms.mjs"))) {
  console.error("usage: node run.mjs <directory of openms.mjs> [--benchmark]");
  process.exit(2);
}
const threadedBuild = JSON.parse(fs.readFileSync(path.join(packageDir, "package.json"), "utf8")).name.endsWith("-mt");

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
  { isolated: false, host: "worker", numThreads: 4 },
  { isolated: false, host: "window" },
  { isolated: true, host: "worker" },
  { isolated: true, host: "worker", numThreads: 1 },
  { isolated: true, host: "window", numThreads: 2 },
];
if (threadedBuild) scenarios.push({ isolated: true, host: "worker", numThreads: 4 });
if (benchmark) {
  for (const n of threadedBuild ? [1, 2, 4] : [1]) {
    scenarios.push({ isolated: true, host: "worker", numThreads: n, benchmark: true });
  }
}

/// The threads a scenario has to report.
function expectedThreads(s, cores) {
  if (!threadedBuild) return { threaded: false, available: false, numThreads: 1 };
  if (!s.isolated) return { threaded: true, available: false, numThreads: 1 };
  return { threaded: true, available: true, numThreads: s.numThreads ?? Math.min(4, Math.ceil(cores / 2)) };
}

//------------------------------------------------------------------------------

const require = createRequire(import.meta.url);
const { chromium } = require(process.env.PLAYWRIGHT_MODULE ?? "playwright");
const browser = await chromium.launch({ executablePath: process.env.CHROMIUM_PATH || undefined });
const servers = { true: await serve(true), false: await serve(false) };
console.log(`${threadedBuild ? "@openms/openms-wasm-mt" : "@openms/openms-wasm"} in ${browser.browserType().name()} ${browser.version()}`);

let failures = 0;
let reference;
for (const s of scenarios) {
  const name = `${s.isolated ? "isolated" : "not isolated"}, ${s.host}, numThreads ${s.numThreads ?? "default"}${s.benchmark ? ", benchmark" : ""}`;
  const problems = [];
  const context = await browser.newContext();
  const page = await context.newPage();
  page.on("pageerror", (e) => problems.push(`page error: ${e.message}`));
  const query = new URLSearchParams({ module: "/pkg/openms.mjs", host: s.host });
  if (s.numThreads) query.set("numThreads", s.numThreads);
  if (s.benchmark) query.set("benchmark", "1");
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
    if (r.featureTable.mz.length !== 13) problems.push(`${r.featureTable.mz.length} features, expected 13`);
    reference ??= r.featureTable;
    if (JSON.stringify(r.featureTable) !== JSON.stringify(reference)) problems.push("features differ from the first scenario");
    const peptides = r.search.hits.map((h) => h.sequence);
    if (r.search.status !== "EXECUTION_OK" || JSON.stringify(peptides) !== JSON.stringify(expectedPeptides)) {
      problems.push(`search: ${r.search.status}, ${JSON.stringify(peptides)}`);
    }
    // in a Web Worker, OpenMS must not block the page (the window host blocks it by design)
    if (s.host === "worker" && outcome.maxStallMs > 500) problems.push(`page blocked for ${outcome.maxStallMs} ms`);
  }

  failures += problems.length ? 1 : 0;
  const t = r?.timings ?? {};
  const details = r ? `threads ${r.info.numThreads}; start ${t.start} ms, features ${t.featureFinding} ms, search ${t.search} ms` +
    (r.benchmark ? `, benchmark search ${t.benchmarkSearch} ms (${r.benchmark.psms} PSMs)` : "") +
    `; page blocked at most ${outcome.maxStallMs} ms` : "";
  console.log(`${problems.length ? "FAIL" : "ok  "} ${name}: ${details}`);
  for (const p of problems) console.log(`       ${p}`);
}

await browser.close();
for (const server of Object.values(servers)) server.close();
console.log(failures ? `${failures} of ${scenarios.length} scenarios failed` : `all ${scenarios.length} scenarios passed`);
process.exit(failures ? 1 : 0);
