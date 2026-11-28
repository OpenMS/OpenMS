// The TypeScript declarations of a build (test openms_wasm_typescript):
//
//   node tests/check_dts.mjs <openms.d.ts> [<openms.d.ts of the other build> [--optional]]
//
// - The declarations of the bindings contain no bigint. Embind passes a long or unsigned long (size_t) as BigInt
//   in the wasm64 build; the bindings take and return them as numbers (bindings/binding_utils.h).
// - With a second file, the openms.d.ts of the wasm32 build for a wasm64 one (or the other way round): the
//   declarations of the bindings are the same, i.e. both builds have the same JavaScript API. --optional: a
//   missing second file is reported, not an error (the other build may not have been built).
//
// Both leave out the interface WasmModule: the raw exports of the module (malloc, operator new, ...), not part of
// the API, whose pointers are BigInt in wasm64.

import fs from "node:fs";

const [file, other, flag] = process.argv.slice(2);
if (!file) {
  console.error("usage: node check_dts.mjs <openms.d.ts> [<openms.d.ts of the other build> [--optional]]");
  process.exit(2);
}

/// The declarations of the bindings: the file without the interface WasmModule
function bindings(path) {
  const text = fs.readFileSync(path, "utf8");
  const withoutExports = text.replace(/^interface WasmModule \{\n[\s\S]*?^\}\n/m, "");
  if (withoutExports === text) throw new Error(`${path}: no interface WasmModule (has Emscripten changed the declarations?)`);
  return withoutExports.split("\n");
}

let failed = false;
const lines = bindings(file);
const bigints = lines.map((line, i) => [i, line]).filter(([, line]) => /\bbigint\b|\bBigInt\b/.test(line));
if (bigints.length) {
  failed = true;
  console.error(`${file}: ${bigints.length} declaration(s) with BigInt (take or return a size_t as number: OMS_FN, OMS_METHOD, OMS_CONSTRUCTOR):`);
  for (const [, line] of bigints) console.error(`  ${line.trim()}`);
} else {
  console.log(`${file}: no BigInt in the declarations of the bindings`);
}

if (other) {
  if (!fs.existsSync(other)) {
    if (flag !== "--optional") {
      console.error(`${other}: not found`);
      process.exit(1);
    }
    console.log(`${other}: not found, not compared (build the other preset to compare the declarations of both builds)`);
  } else {
    const otherLines = bindings(other);
    const n = Math.max(lines.length, otherLines.length);
    const differences = [];
    for (let i = 0; i < n; ++i) {
      if (lines[i] !== otherLines[i]) differences.push(i);
    }
    if (differences.length) {
      failed = true;
      console.error(`${file} and ${other} declare different bindings (${differences.length} line(s)), e.g.:`);
      for (const i of differences.slice(0, 10)) {
        console.error(`  line ${i + 1}:\n    ${lines[i] ?? "(none)"}\n    ${otherLines[i] ?? "(none)"}`);
      }
    } else {
      console.log(`${other}: the same declarations of the bindings`);
    }
  }
}
process.exit(failed ? 1 : 0);
