// The loader of the assembled package @openms/openms-wasm (loader/index.mjs, assemble_package.cmake):
// it starts the threaded build where threads are available (Node.js: always), else the single-threaded one.
//
// OPENMS_WASM_PACKAGE  directory of the assembled package (the tests are skipped without it)

import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { pathToFileURL } from "node:url";

const packageDir = process.env.OPENMS_WASM_PACKAGE;
const skip = !packageDir && "OPENMS_WASM_PACKAGE is not set";
const quiet = { print: () => {}, printErr: () => {} };

async function loader() {
  return import(pathToFileURL(path.join(packageDir, "index.mjs")).href);
}

test("package layout", { skip }, () => {
  for (const file of ["index.mjs", "index.d.ts", "openms.d.ts", "package.json", "st/openms.mjs", "st/openms.wasm",
                      "mt/openms.mjs", "mt/openms.wasm"]) {
    assert.ok(fs.existsSync(path.join(packageDir, file)), file);
  }
  const pkg = JSON.parse(fs.readFileSync(path.join(packageDir, "package.json"), "utf8"));
  assert.equal(pkg.name, "@openms/openms-wasm");
  assert.equal(pkg.exports["."].default, "./index.mjs");
});

test("Node.js gets the threaded build", { skip }, async () => {
  const { default: createOpenMS, threadSupport } = await loader();
  assert.deepEqual(threadSupport(), { supported: true, reason: "Node.js worker threads" });
  const m = await createOpenMS({ ...quiet, numThreads: 2 });
  assert.deepEqual(m.getThreadingInfo(), { threaded: true, available: true, numThreads: 2 });
});

test("threads: false gives the single-threaded build", { skip }, async () => {
  const { default: createOpenMS } = await loader();
  const m = await createOpenMS({ ...quiet, threads: false, numThreads: 4 });
  assert.deepEqual(m.getThreadingInfo(), { threaded: false, available: false, numThreads: 1 });
  // the same module: a quick computation
  const seq = m.AASequence.fromString("PEPTIDE");
  assert.ok(Math.abs(seq.getMonoWeight() - 799.35996) < 1e-4);
  seq.delete();
});

test("invalid threads option", { skip }, async () => {
  const { default: createOpenMS } = await loader();
  await assert.rejects(createOpenMS({ threads: "yes" }), /threads must be "auto", true or false/);
});
