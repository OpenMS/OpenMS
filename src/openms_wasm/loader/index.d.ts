// TypeScript declarations of the entry point of @openms/openms-wasm (index.mjs). The module
// itself, the same for both builds, is declared in openms.d.ts.

import type { MainModule } from "./openms";

export * from "./openms";

export interface OpenMSOptions {
  /** "auto" (default): threads where the environment allows them; true / false: that build. */
  threads?: boolean | "auto";
  /** Number of threads of the threaded build (default: Node.js: OMP_NUM_THREADS or all cores; browsers: half the cores, at most 4). */
  numThreads?: number;
  /** Further options of the Emscripten module (print, printErr, locateFile, ...). */
  [option: string]: unknown;
}

export interface ThreadSupport {
  /** Whether the threaded build can run with more than one thread here. */
  supported: boolean;
  /** Why (not): e.g. that the page is not cross-origin isolated. */
  reason: string;
}

/** Whether threads are available here, and if not, why. */
export function threadSupport(): ThreadSupport;

/** Starts OpenMS: the threaded build where threads are available, else the single-threaded one. */
export default function createOpenMS(options?: OpenMSOptions): Promise<MainModule>;
