// Reading and writing files: through the virtual file system and from strings.

import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { pathToFileURL } from "node:url";
import { loadOpenMS, using, assertClose, testDataDir, modulePath } from "./helpers.mjs";

const OpenMS = await loadOpenMS();

/// MzMLFile_test: MzMLFile_1.mzML has 4 spectra (15, 10, 15 and 0 peaks) and 2 chromatograms
function checkMzMLFile1(exp) {
  assert.equal(exp.getNrSpectra(), 4);
  assert.equal(exp.getNrChromatograms(), 2);
  using([exp.getSpectrum(0), exp.getSpectrum(1)], (s0, s1) => {
    assert.equal(s0.size(), 15);
    assert.equal(s1.size(), 10);
    assert.equal(s0.getMSLevel(), 1);
    assert.equal(s1.getMSLevel(), 2);
    assertClose(assert, s0.getRT(), 5.1, 1e-6, "RT");
    const { mz, intensity } = s0.getPeaks();
    assert.ok(mz instanceof Float64Array && intensity instanceof Float32Array);
    assertClose(assert, mz[0], 0, 1e-6, "first m/z");
    assertClose(assert, intensity[0], 15, 1e-6, "first intensity");
    const precursors = s1.getPrecursors();
    assert.equal(precursors.length, 2);
    // the selected ion m/z (5.5555); MzMLFile_test compares it to the isolation window target 5.55 with tolerance 0.01
    assertClose(assert, precursors[0].mz, 5.5555, 1e-6, "precursor m/z");
    assert.equal(precursors[0].charge, 2);
    assertClose(assert, precursors[0].intensity, 120053, 1e-3, "precursor intensity");
    assertClose(assert, precursors[0].driftTime, 8.1, 1e-6, "precursor drift time");
    assert.ok(precursors[0].activationMethods.includes("CID"));
    assert.equal(s0.getNativeID(), "index=0");
    assert.equal(s0.getMetaValue("maldi_spot_id"), "M0");
    assert.equal(s0.getType(), "Centroid");
  });
}

test("mountHostDirectory and unmount", () => {
  OpenMS.mountHostDirectory("/mounted", testDataDir);
  assert.ok(OpenMS.FS.stat("/mounted/MzMLFile_1.mzML").size > 1000);
  assert.ok(OpenMS.FS.readdir("/mounted").includes("MzMLFile_1.mzML"));
  OpenMS.unmount("/mounted");
  assert.throws(() => OpenMS.FS.stat("/mounted/MzMLFile_1.mzML"));
  assert.throws(() => OpenMS.unmount("/mounted"), /mounted/);
  assert.throws(() => OpenMS.mountHostDirectory("/no/such/parent/x", testDataDir), /no\/such\/parent/);
});

/// Descriptor @p fd of the process, opened on @p file, if the lowest free descriptors lead there; else undefined
/// (the descriptors taken on the way are closed again).
function holdDescriptor(fd, file) {
  const taken = [];
  try {
    for (;;) {
      const next = fs.openSync(file, "r");
      if (next === fd) return next;
      taken.push(next);
      if (next > fd) return undefined;
    }
  } finally {
    for (const t of taken) fs.closeSync(t);
  }
}

/// Whether descriptor @p fd of the process is open.
function isOpen(fd) {
  try {
    fs.fstatSync(fd);
    return true;
  } catch {
    return false;
  }
}

// WasmFS' Node.js backend took the errno of a failed open of a host file for a file descriptor of the process
// (wasmfs_node.js): a file deleted on the host after the module had looked it up was open on descriptor 44
// (ENOENT), which the module read and then closed, whatever it was.
test("mountHostDirectory: a file deleted on the host cannot be read, and no descriptor of the process is closed", () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "openms-wasm-io-"));
  const file = path.join(dir, "deleted.txt");
  fs.writeFileSync(file, "deleted on the host\n");
  // descriptor 44 of the test's own, where the lowest free descriptors lead there (in the threaded module, the
  // Workers of its threads may have it)
  const own = holdDescriptor(44, file);
  let readError;
  let ownStillOpen;
  OpenMS.mountHostDirectory("/host", dir);
  try {
    assert.equal(OpenMS.FS.readFile("/host/deleted.txt", { encoding: "utf8" }), "deleted on the host\n");
    fs.unlinkSync(file);
    try {
      OpenMS.FS.readFile("/host/deleted.txt");
    } catch (e) {
      readError = e;
    }
    if (own !== undefined) ownStillOpen = isOpen(own);
  } finally {
    if (own !== undefined && isOpen(own)) fs.closeSync(own);
    OpenMS.unmount("/host");
    fs.rmSync(dir, { recursive: true, force: true });
  }
  if (own !== undefined) assert.ok(ownStillOpen, "descriptor 44 of the test is still open");
  assert.equal(readError?.errno, 44, "the read fails with ENOENT");
});

/// The open descriptors of the process, number -> what they point to (Linux: /proc/self/fd).
function descriptorTable() {
  const table = new Map();
  for (const name of fs.readdirSync("/proc/self/fd")) {
    try {
      table.set(Number(name), fs.readlinkSync(`/proc/self/fd/${name}`));
    } catch {
      // the descriptor of the readdir itself, closed again
    }
  }
  return table;
}

// The other half of the same bug (wasmfs_node.js): a close of a host file that fails (here EBADF, 8, because the
// descriptor was closed behind the module's back) was reported as a positive number, which came out of
// FS.close() as the errno 65528.
//
// Closing the module's descriptor on the host frees its number for a moment, and the module closes that number
// again: the test therefore has to stay before any test that starts threads of the module (their Workers open and
// close descriptors of the process, one of them could take the number). It checks that the number is free right
// before FS.close(), and that no other descriptor of the process changed.
test("mountHostDirectory: a host file that cannot be closed reports its errno", { skip: !fs.existsSync("/proc/self/fd") }, () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "openms-wasm-io-"));
  const file = path.join(dir, "closed.txt");
  fs.writeFileSync(file, "closed behind the module's back\n");
  OpenMS.mountHostDirectory("/host", dir);
  try {
    const stream = OpenMS.FS.open("/host/closed.txt", "r");
    const before = descriptorTable();
    const target = fs.realpathSync(file);
    const fd = [...before].find(([, pointsTo]) => pointsTo === target)?.[0];
    assert.notEqual(fd, undefined, "the module has the file open");
    fs.closeSync(fd);
    const expected = new Map(before);
    expected.delete(fd);
    assert.ok(!isOpen(fd), `descriptor ${fd} was taken again before the module closed it`);
    assert.throws(() => OpenMS.FS.close(stream), (e) => e.errno === 8, "FS.close() fails with EBADF");
    assert.deepEqual(descriptorTable(), expected, "no other descriptor of the process was closed or opened");
  } finally {
    OpenMS.unmount("/host");
    fs.rmSync(dir, { recursive: true, force: true });
  }
});

test("writeFileChunk, readFileChunk and allocateFile: files in parts, byte for byte", () => {
  const bytes = new Uint8Array(5 * 1048576 + 12345);
  for (let i = 0; i < bytes.length; ++i) bytes[i] = (i * 2654435761) >>> 24;
  const chunk = 1048576;
  OpenMS.FS.mkdir("/chunks");
  try {
    // allocated first: written in place
    OpenMS.allocateFile("/chunks/a.bin", bytes.length);
    assert.equal(OpenMS.FS.stat("/chunks/a.bin").size, bytes.length);
    for (let at = 0; at < bytes.length; at += chunk) {
      OpenMS.writeFileChunk("/chunks/a.bin", bytes.subarray(at, at + chunk), at);
    }
    assert.deepEqual(OpenMS.FS.readFile("/chunks/a.bin"), bytes);
    // not allocated: created, and it grows; chunks in any order
    const order = [3, 0, 5, 1, 4, 2];
    for (const k of order) OpenMS.writeFileChunk("/chunks/b.bin", bytes.subarray(k * chunk, (k + 1) * chunk), k * chunk);
    assert.deepEqual(OpenMS.FS.readFile("/chunks/b.bin"), bytes);
    // read back in parts: fewer bytes at the end, none past it
    const parts = [];
    for (let at = 0; at < bytes.length; at += chunk) {
      const part = OpenMS.readFileChunk("/chunks/a.bin", at, chunk);
      assert.ok(part instanceof Uint8Array);
      parts.push(part);
    }
    assert.equal(parts.at(-1).length, bytes.length % chunk);
    assert.deepEqual(new Uint8Array(Buffer.concat(parts)), bytes);
    assert.equal(OpenMS.readFileChunk("/chunks/a.bin", bytes.length, 10).length, 0);
    assert.equal(OpenMS.readFileChunk("/chunks/a.bin", bytes.length + 100, 10).length, 0);
    assert.deepEqual(OpenMS.readFileChunk("/chunks/a.bin", 7, 3), bytes.subarray(7, 10));
    // FS.write() (one copy instead of byte by byte, wasmfs_fs.js) with an offset into its buffer
    const stream = OpenMS.FS.open("/chunks/c.bin", "w");
    OpenMS.FS.write(stream, bytes, 100, 1000, 0);
    OpenMS.FS.close(stream);
    assert.deepEqual(OpenMS.FS.readFile("/chunks/c.bin"), bytes.subarray(100, 1100));
    // allocateFile resizes, keeping the contents up to the new size
    OpenMS.allocateFile("/chunks/a.bin", 10);
    assert.deepEqual(OpenMS.FS.readFile("/chunks/a.bin"), bytes.subarray(0, 10));
    // errors
    assert.throws(() => OpenMS.writeFileChunk("/chunks/a.bin", [1, 2, 3], 0), /Uint8Array/);
    assert.throws(() => OpenMS.writeFileChunk("/chunks/a.bin", bytes.subarray(0, 1), -1), /offset/);
    assert.throws(() => OpenMS.writeFileChunk("/chunks/a.bin", bytes.subarray(0, 1), 0.5), /offset/);
    assert.throws(() => OpenMS.readFileChunk("/chunks/missing.bin", 0, 1), /FileNotFound.*missing\.bin/);
    assert.throws(() => OpenMS.writeFileChunk("/no/such/dir/x.bin", bytes.subarray(0, 1), 0), /UnableToCreateFile.*x\.bin/);
    assert.throws(() => OpenMS.allocateFile("/chunks/a.bin", -5), /size/);
  } finally {
    for (const f of ["a.bin", "b.bin", "c.bin"]) {
      try { OpenMS.FS.unlink(`/chunks/${f}`); } catch { /* not there */ }
    }
    OpenMS.FS.rmdir("/chunks");
  }
});

test("allocateFile and writeFileChunk: a file has to be smaller than 2 GiB", () => {
  // WasmFS' memory backend keeps a file in one std::vector: less than 2 GiB in wasm32, and a block close to that cannot
  // be allocated. Sizes beyond the limit are rejected before any allocation, also where they exceed what a size_t of
  // wasm32 holds (4 GiB).
  const MiB = 1048576, limit = 2 ** 31 - 64 * MiB;
  const tooLarge = /InvalidValue: the value '\d+' was used but is not valid; a file of the module's file system has to be smaller than 2 GiB \(at most 2080374784 bytes\): \/limit\.bin/;
  try {
    for (const size of [limit + 1, 2 ** 31 - 1, 2 ** 31, 2 ** 32 + 5, 2 ** 40]) {
      assert.throws(() => OpenMS.allocateFile("/limit.bin", size), tooLarge, `allocateFile ${size}`);
    }
    const data = new Uint8Array(16);
    for (const offset of [limit - 8, 2 ** 31, 2 ** 32, 2 ** 32 - 8]) {
      assert.throws(() => OpenMS.writeFileChunk("/limit.bin", data, offset), tooLarge, `writeFileChunk at ${offset}`);
    }
    // up to the limit: a small file is fine
    OpenMS.writeFileChunk("/limit.bin", data, 0);
    assert.equal(OpenMS.FS.stat("/limit.bin").size, 16);
  } finally {
    try { OpenMS.FS.unlink("/limit.bin"); } catch { /* not there */ }
  }
});

test("FS.writeFile, FS.write, FS.read and the chunk functions fail cleanly where the memory cannot hold their data", async () => {
  // Emscripten's FS.writeFile/write/read took the 0 that malloc returns then for the address of their buffer, and
  // wrote over the start of the module's memory (wasmfs_fs.js). A module of its own: its memory grows to the limit.
  const { default: createOpenMS } = await import(pathToFileURL(modulePath).href);
  const m = await createOpenMS({ numThreads: 1, print: () => {}, printErr: () => {} });
  const MB = 1048576;
  const { heapMax } = m.getMemoryInfo();
  // files that take 3/4 of the memory: one of the largest size the memory file system takes (2 GiB - 64 MiB; one byte
  // more is rejected) and one of 1 GiB; then data of 1/3 of the memory, which the memory has (it is larger), but not free
  const occupying = ["/max.bin", "/occupied.bin"];
  m.allocateFile(occupying[0], 2 ** 31 - 64 * MB);
  assert.equal(m.FS.stat(occupying[0]).size, 2 ** 31 - 64 * MB);
  assert.throws(() => m.allocateFile(occupying[0], 2 ** 31 - 64 * MB + 1), /smaller than 2 GiB/);
  m.allocateFile(occupying[1], 1024 * MB);
  const tooLarge = Math.floor(heapMax / 3);
  const data = new Uint8Array(tooLarge); // (never written: its pages stay unallocated)
  const enomem = (e) => e?.name === "ErrnoError" && e.errno === 48;
  assert.throws(() => m.FS.writeFile("/x.bin", data), enomem, "FS.writeFile");
  const stream = m.FS.open("/y.bin", "w+");
  assert.throws(() => m.FS.write(stream, data, 0, data.length, 0), enomem, "FS.write");
  assert.throws(() => m.FS.read(stream, data, 0, data.length, 0), enomem, "FS.read");
  m.FS.close(stream);
  assert.throws(() => m.allocateFile("/z.bin", tooLarge), /OutOfMemory/, "allocateFile");
  assert.throws(() => m.writeFileChunk("/z.bin", data, 0), /OutOfMemory/, "writeFileChunk");
  // growing a file of 1 GiB takes a block of 2 GiB (it grows by doubling), which the file system cannot take
  assert.throws(() => m.writeFileChunk(occupying[1], data.subarray(0, 16), 1024 * MB),
    /growing the file \/occupied\.bin of 1073741824 bytes takes a block of 2147483647 bytes.*allocate the file at its final size first/);
  // FS.write() and FS.read() of a range beyond their buffer: EINVAL, before any allocation
  const einval = (e) => e?.name === "ErrnoError" && e.errno === 28;
  const small = m.FS.open("/small.bin", "w+");
  assert.throws(() => m.FS.write(small, new Uint8Array(4), 2, 4, 0), einval, "FS.write past its buffer");
  assert.throws(() => m.FS.read(small, new Uint8Array(4), 0, 8, 0), einval, "FS.read past its buffer");
  m.FS.write(small, new Uint8Array([1, 2, 3, 4]), 1, 2, 0);
  m.FS.close(small);
  assert.deepEqual(m.FS.readFile("/small.bin"), new Uint8Array([2, 3]));
  // growing a file by a chunk: the memory has to hold it at its new size (here: twice its size)
  m.allocateFile("/w.bin", Math.floor(tooLarge / 2));
  assert.throws(() => m.writeFileChunk("/w.bin", data.subarray(0, 16), Math.floor(tooLarge / 2)), /OutOfMemory/, "writeFileChunk past the end");
  // the module is usable: the start of its memory is intact, and the file system works
  for (const file of occupying) m.FS.unlink(file);
  m.FS.unlink("/w.bin");
  m.FS.writeFile("/ok.bin", new Uint8Array([1, 2, 3]));
  assert.deepEqual(m.FS.readFile("/ok.bin"), new Uint8Array([1, 2, 3]));
  assert.equal(m.FS.stat("/y.bin").size, 0);
  assert.match(m.getVersion(), /^\d+\.\d+/);
  using([m.AASequence.fromString("PEPTIDEK")], (seq) => assertClose(assert, seq.getMonoWeight(m.ResidueType.Full, 0), 927.4549, 1e-3, "mono weight"));
});

test("MzMLFile from the (mounted) file system", () => {
  using([new OpenMS.MzMLFile(), new OpenMS.MSExperiment()], (file, exp) => {
    file.load("/data/MzMLFile_1.mzML", exp);
    checkMzMLFile1(exp);
    const summary = exp.getSpectrumSummary();
    assert.equal(summary.length, 4);
    assert.equal(summary[1].msLevel, 2);
    assert.deepEqual(exp.getMSLevels(), [1, 2]);
  });
});

test("MzMLFile from a string, and round trip through MEMFS", () => {
  const text = fs.readFileSync(path.join(testDataDir, "MzMLFile_1.mzML"), "utf8");
  using([new OpenMS.MzMLFile(), new OpenMS.MSExperiment(), new OpenMS.MSExperiment(), new OpenMS.MSExperiment()],
    (file, exp, fromFile, fromString) => {
      file.loadBuffer(text, exp);
      checkMzMLFile1(exp);

      file.store("/tmp/roundtrip.mzML", exp);
      assert.ok(OpenMS.FS.stat("/tmp/roundtrip.mzML").size > 1000);
      file.load("/tmp/roundtrip.mzML", fromFile);
      checkMzMLFile1(fromFile);

      file.loadBuffer(file.storeBuffer(exp), fromString);
      checkMzMLFile1(fromString);
    });
});

test("MzMLFile load options", () => {
  using([new OpenMS.MzMLFile(), new OpenMS.MSExperiment()], (file, exp) => {
    file.setOptions({ msLevels: [2] });
    assert.deepEqual(file.getOptions().msLevels, [2]);
    file.load("/data/MzMLFile_1.mzML", exp);
    assert.deepEqual(exp.getMSLevels(), [2]);
  });
});

test("FileHandler", () => {
  assert.equal(OpenMS.FileHandler.getType("/data/MzMLFile_1.mzML"), "mzML");
  assert.equal(OpenMS.FileHandler.getTypeByFileName("x.featureXML"), "featureXML");
  using([new OpenMS.FileHandler(), new OpenMS.MSExperiment(), new OpenMS.FeatureMap()], (handler, exp, features) => {
    handler.loadExperiment("/data/MzMLFile_1.mzML", exp);
    checkMzMLFile1(exp);
    handler.loadFeatures("/data/FeatureXMLFile_1.featureXML", features);
    assert.equal(features.size(), 2);
  });
});

test("FASTAFile", () => {
  using([new OpenMS.FASTAFile()], (file) => {
    const entries = file.load("/data/FASTAFile_test.fasta");
    assert.ok(entries.length > 0);
    assert.ok(entries.every((e) => typeof e.identifier === "string" && e.sequence.length > 0));
    file.store("/tmp/out.fasta", entries);
    assert.deepEqual(file.load("/tmp/out.fasta"), entries);
  });
});

test("IdXMLFile", () => {
  using([new OpenMS.IdXMLFile(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (file, proteins, peptides) => {
      file.load("/data/IdXMLFile_whole.idXML", proteins, peptides);
      // IdXMLFile_test: 2 protein and 3 peptide identification runs
      assert.equal(proteins.size(), 2);
      assert.equal(peptides.size(), 3);
      const ids = peptides.toArray();
      assert.equal(ids[0].hits.length, 2);
      assert.equal(ids[0].hits[0].sequence, "PEPTIDER");
      assertClose(assert, ids[0].hits[0].score, 0.9, 1e-9, "score");
      assert.equal(ids[0].hits[0].charge, 1);
      const runs = proteins.toArray();
      assert.equal(runs[0].searchEngine, "Mascot");
      assert.equal(runs[0].scoreType, "MOWSE");
      assert.equal(runs[0].metaValues.name, "ProteinIdentification");
      using([peptides.get(0)], (id) => {
        assert.equal(id.getHitCount(), 2);
        const hits = id.getHits();
        try {
          assert.equal(hits[0].getSequenceString(), "PEPTIDER");
        } finally {
          hits.forEach((h) => h.delete());
        }
      });
      file.store("/tmp/out.idXML", proteins, peptides);
      using([new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()], (p2, pep2) => {
        file.load("/tmp/out.idXML", p2, pep2);
        assert.deepEqual(pep2.toArray().map((i) => i.hits.map((h) => h.sequence)),
                         ids.map((i) => i.hits.map((h) => h.sequence)));
      });
    });
});

test("FeatureXMLFile", () => {
  using([new OpenMS.FeatureXMLFile(), new OpenMS.FeatureMap()], (file, features) => {
    file.load("/data/FeatureXMLFile_1.featureXML", features);
    assert.equal(features.size(), 2);
    const table = features.getFeatureTable();
    assert.equal(table.rt.length, 2);
    using([features.get(0)], (f) => {
      assertClose(assert, f.getRT(), table.rt[0], 1e-9, "RT");
      assertClose(assert, f.getMZ(), table.mz[0], 1e-9, "m/z");
    });
  });
});

test("MSExperiment.getSpectrumSummary: m/z range of peaks not sorted by m/z", () => {
  using([new OpenMS.MSExperiment(), new OpenMS.MSSpectrum()], (exp, spectrum) => {
    spectrum.setPeaks([500, 100, 300], [1, 2, 3]);
    exp.addSpectrum(spectrum);
    const summary = exp.getSpectrumSummary();
    assert.equal(summary.length, 1);
    assert.equal(summary[0].minMZ, 100);
    assert.equal(summary[0].maxMZ, 500);
  });
});

test("MSExperiment.get2DPeakData: MS level 0 is all levels", () => {
  using([new OpenMS.MSExperiment(), new OpenMS.MSSpectrum(), new OpenMS.MSSpectrum()], (exp, ms1, ms2) => {
    ms1.setRT(1);
    ms1.setMSLevel(1);
    ms1.setPeaks([100], [10]);
    ms2.setRT(2);
    ms2.setMSLevel(2);
    ms2.setPeaks([200], [20]);
    exp.addSpectrum(ms1);
    exp.addSpectrum(ms2);
    // (minRT, maxRT, minMZ, maxMZ, msLevel)
    const all = exp.get2DPeakData(0, 10, 0, 1000, 0);
    assert.deepEqual(Array.from(all.rt), [1, 2]);
    assert.deepEqual(Array.from(all.mz), [100, 200]);
    const ms2Only = exp.get2DPeakData(0, 10, 0, 1000, 2);
    assert.deepEqual(Array.from(ms2Only.rt), [2]);
    assert.deepEqual(Array.from(ms2Only.mz), [200]);
  });
});
