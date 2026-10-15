// Thread settings of the multi-threaded OpenMS module (OPENMS_WASM_THREADS), evaluated
// before the runtime starts (linked with --pre-js; Module is the object passed to the
// factory, createOpenMS()).
//
// Module.numThreads: the number of threads OpenMS computes with (OpenMP's OMP_NUM_THREADS).
//  - Node.js: threads start when they are needed. Default: OMP_NUM_THREADS of the process,
//    else all logical cores, like native OpenMS.
//  - Browsers: threads share the module's memory, which a page only may when it is
//    cross-origin isolated (crossOriginIsolated: served with COOP/COEP headers). Elsewhere
//    the module runs single-threaded, whatever numThreads says. The threads are Web Workers,
//    started with the module (openmsThreadPoolSize, Emscripten's PTHREAD_POOL_SIZE): a
//    thread that waits for another one blocks, and a Worker cannot start while its creator
//    blocks, so a team of OpenMP threads cannot wait for new Workers. Default: half the
//    logical cores, at most 4 (as ONNX Runtime Web), leaving the rest to the page.

/** OpenMP team size; undefined: OpenMP's default (all logical cores). */
var openmsNumThreads = (() => {
  const requested = Module["numThreads"];
  if (requested !== undefined && !(Number.isInteger(requested) && requested >= 1)) {
    throw new Error(`numThreads must be a positive integer, got ${requested}`);
  }
  if (typeof process == "object" && process.versions?.node) {
    return requested;
  }
  if (!globalThis.crossOriginIsolated) {
    return 1;
  }
  const cores = globalThis.navigator?.hardwareConcurrency || 1;
  return requested ?? Math.min(4, Math.ceil(cores / 2));
})();

/**
 * Workers started with the module. Node.js: none (it starts them on demand). Browsers: one
 * per thread of a team but the calling one; three for the helper threads OpenMS starts besides
 * OpenMP (e.g. the writers of OpenSwath); and the OpenMP team of such a helper thread: ProSE
 * reads the spectra on a helper thread with a quarter of the threads (at least 2, at most 16)
 * while the calling thread's team builds the index. A thread that finds the pool empty waits
 * for a new Worker, which cannot start while its creator blocks: the search would hang.
 */
var openmsThreadPoolSize = (() => {
  if ((typeof process == "object" && process.versions?.node) || openmsNumThreads == 1) return 0;
  const helperTeam = Math.min(Math.max(Math.floor(openmsNumThreads / 4), Math.min(openmsNumThreads, 2)), 16);
  return (openmsNumThreads - 1) + 3 + (helperTeam - 1);
})();

/**
 * In browsers, the threads' Web Workers are started from a Blob that imports this module, not
 * from the module's URL: a Worker started from a Blob inherits the cross-origin isolation of its
 * creator, while one started from a URL needs the Cross-Origin-Embedder-Policy header on that
 * response. The module usually runs in a Web Worker of the page, and Chrome fetches the script
 * of a Worker created inside a Worker directly from the server, past a service worker that
 * would add the header (examples/sw.js): on a server without the header, the threads would
 * not start. (Node.js workers are started from the module's path as usual.)
 */
if (!(typeof process == "object" && process.versions?.node) && Module["mainScriptUrlOrBlob"] === undefined) {
  Module["mainScriptUrlOrBlob"] = new Blob([`import ${JSON.stringify(import.meta.url)};`], { type: "text/javascript" });
}

Module["preRun"] = [].concat(Module["preRun"] || [], () => {
  // The OpenMP variables of the Node.js process (Emscripten's environment starts empty).
  if (typeof process == "object" && process.versions?.node) {
    for (const [name, value] of Object.entries(process.env)) {
      if (/^(OMP|KMP)_/.test(name)) ENV[name] = value;
    }
  }
  // Process-wide: also the default of threads OpenMS starts itself.
  if (openmsNumThreads !== undefined) ENV["OMP_NUM_THREADS"] = String(openmsNumThreads);
});
