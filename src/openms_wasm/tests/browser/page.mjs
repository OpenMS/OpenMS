// Runs the workload on this page or in a Web Worker and publishes the outcome as
// window.testResult (read by run.mjs). While the workload runs, a 10 ms timer measures how
// long the page cannot react (maxStallMs).
import { runWorkload } from "./workload.mjs";

const params = new URLSearchParams(location.search);
const job = {
  moduleUrl: new URL(params.get("module"), location.href).href,
  dataUrl: new URL("/data", location.href).href,
  numThreads: params.has("numThreads") ? Number(params.get("numThreads")) : undefined,
  benchmark: params.get("benchmark") === "1",
  fileByFile: params.get("fileByFile") === "1",
  manyFiles: Number(params.get("manyFiles") ?? 0),
};
const page = { crossOriginIsolated: globalThis.crossOriginIsolated, sharedArrayBuffer: typeof SharedArrayBuffer != "undefined",
               hardwareConcurrency: navigator.hardwareConcurrency };

let last = performance.now(), maxStallMs = 0;
const heartbeat = setInterval(() => {
  const now = performance.now();
  maxStallMs = Math.max(maxStallMs, now - last - 10);
  last = now;
}, 10);

function finish(outcome) {
  clearInterval(heartbeat);
  window.testResult = { page, host: params.get("host"), maxStallMs: Math.round(maxStallMs), ...outcome };
}

if (params.get("host") === "worker") {
  const worker = new Worker(new URL("worker.mjs", import.meta.url), { type: "module" });
  worker.onmessage = ({ data }) => finish(data);
  worker.onerror = (e) => finish({ error: `worker error: ${e.message}` });
  worker.postMessage(job);
} else {
  // give the timer a first tick before the page blocks
  setTimeout(() => runWorkload(job).then((result) => finish({ result }), (e) => finish({ error: String(e?.stack ?? e) })), 20);
}
