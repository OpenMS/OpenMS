// OpenMS for JavaScript: the entry point of @openms/openms-wasm. The package holds two
// builds of the module with the same API, st/ (single-threaded) and mt/ (with threads:
// OpenMP and the helper threads of OpenMS on Web Workers or Node.js worker threads), and
// this loader starts the one the environment can run.
//
//   import createOpenMS from "@openms/openms-wasm";
//   const OpenMS = await createOpenMS();
//   OpenMS.getThreadingInfo();   // { threaded, available, numThreads }

const isNode = typeof process == "object" && !!process.versions?.node;

/**
 * Whether the threaded build can run with more than one thread here, and if not, why.
 * Threads share the module's memory; a browser allows that only on a cross-origin isolated
 * page, served with the headers Cross-Origin-Opener-Policy: same-origin and
 * Cross-Origin-Embedder-Policy: require-corp (or credentialless: Chrome, Firefox).
 */
export function threadSupport() {
  if (isNode) {
    return { supported: true, reason: "Node.js worker threads" };
  }
  if (globalThis.crossOriginIsolated !== true) {
    return {
      supported: false,
      reason: "the page is not cross-origin isolated (serve it with the headers " +
              "Cross-Origin-Opener-Policy: same-origin and Cross-Origin-Embedder-Policy: require-corp)",
    };
  }
  if (typeof Worker == "undefined") {
    return { supported: false, reason: "this context cannot start Web Workers (Safari before 16.4 in a Web Worker)" };
  }
  try {
    // shared memory has to reach the threads (Firefox can isolate a page and still refuse this)
    new MessageChannel().port1.postMessage(new WebAssembly.Memory({ initial: 1, maximum: 1, shared: true }));
  } catch (e) {
    return { supported: false, reason: `shared WebAssembly memory cannot be passed to workers (${e?.message ?? e})` };
  }
  return { supported: true, reason: "the page is cross-origin isolated" };
}

/**
 * Starts OpenMS.
 * @param options.threads    "auto" (default): the threaded build where threadSupport() allows it,
 *                           else the single-threaded one; true / false: that build. (The threaded
 *                           build also runs where threads are unavailable, on one thread.)
 * @param options.numThreads the number of threads of the threaded build (default: Node.js:
 *                           OMP_NUM_THREADS, else all cores; browsers: half the cores, at most 4)
 * Other options go to the module (print, printErr, locateFile, ...).
 */
export default async function createOpenMS({ threads = "auto", ...options } = {}) {
  if (threads !== "auto" && typeof threads != "boolean") {
    throw new Error(`threads must be "auto", true or false, got ${threads}`);
  }
  const support = threadSupport();
  if (threads !== false && !support.supported && options.numThreads > 1) {
    console.warn(`OpenMS: numThreads is ${options.numThreads}, but OpenMS runs on one thread: ${support.reason}.`);
  }
  // Two literal imports, so that bundlers keep both builds (only the one used is downloaded).
  const { default: factory } = threads === true || (threads === "auto" && support.supported)
    ? await import("./mt/openms.mjs")
    : await import("./st/openms.mjs");
  return factory(options);
}
