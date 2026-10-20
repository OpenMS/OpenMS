// Shared setup of the tests of the OpenMS WebAssembly module (run: node --test tests/).
//
// OPENMS_WASM_MODULE  path of openms.mjs (default: the build tree of the preset wasm32-emscripten-release)
// OPENMS_TEST_DATA    the class test data directory (default: src/tests/class_tests/openms/data)

import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const sourceRoot = path.resolve(here, "../../..");

export const modulePath = process.env.OPENMS_WASM_MODULE
  ?? path.join(sourceRoot, "build/wasm32-emscripten-release/openms-wasm/openms.mjs");
export const testDataDir = process.env.OPENMS_TEST_DATA
  ?? path.join(sourceRoot, "src/tests/class_tests/openms/data");
export const toppDataDir = path.resolve(testDataDir, "../../../topp");

let instance;

/// The module, instantiated once per test file; the test data directories are mounted at
/// /data and /topp of its virtual file system.
export async function loadOpenMS() {
  if (!instance) {
    const { default: createOpenMS } = await import(pathToFileURL(modulePath).href);
    instance = await createOpenMS({ print: () => {}, printErr: () => {} });
    instance.mountHostDirectory("/data", testDataDir);
    instance.mountHostDirectory("/topp", toppDataDir);
  }
  return instance;
}

/// Runs @p fn with objects that are deleted afterwards (Embind objects are not garbage collected).
export function using(objects, fn) {
  try {
    return fn(...objects);
  } finally {
    for (const o of objects) o.delete();
  }
}

/// |a - b| <= tolerance
export function assertClose(assert, actual, expected, tolerance, message) {
  assert.ok(Math.abs(actual - expected) <= tolerance,
    `${message ?? "value"}: expected ${expected} ± ${tolerance}, got ${actual}`);
}
