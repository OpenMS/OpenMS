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
 * per thread of a team but the calling one, plus three for the helper threads OpenMS starts
 * besides OpenMP (e.g. the writers of OpenSwath).
 */
var openmsThreadPoolSize = (typeof process == "object" && process.versions?.node) || openmsNumThreads == 1
  ? 0 : openmsNumThreads + 2;

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
