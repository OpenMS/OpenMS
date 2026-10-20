# OpenMS for WebAssembly

The OpenMS C++ library compiled to WebAssembly, with JavaScript bindings. Web
applications (and Node.js programs) can use OpenMS directly, without native binaries,
a server or a Python runtime: the library runs in the browser, on the user's data,
which never leaves the machine.

```js
import createOpenMS from "@openms/openms-wasm";

const OpenMS = await createOpenMS();
const peptide = OpenMS.AASequence.fromString("PEPTM(Oxidation)IDEK");
console.log(peptide.getMonoWeight(), peptide.getMZ(2));
peptide.delete();
```

The bindings cover a selection of the OpenMS API, named as in C++ and pyOpenMS:

| Area | Classes and functions |
|---|---|
| Spectra and maps | `MSSpectrum`, `MSChromatogram`, `MSExperiment` (peaks as typed arrays, `getSpectrumSummary()`, `get2DPeakData()`, `calculateTIC()`) |
| Files | `MzMLFile` (also from and to strings), `MzXMLFile`, `MascotGenericFile`, `FASTAFile`, `FeatureXMLFile`, `ConsensusXMLFile`, `IdXMLFile`, `MzIdentMLFile`, `FileHandler` |
| Chemistry | `AASequence`, `EmpiricalFormula`, `IsotopeDistribution`, `CoarseIsotopePatternGenerator`, `FineIsotopePatternGenerator`, `ResidueDB`, `ModificationsDB`, `ProteaseDB`, `ProteaseDigestion`, `TheoreticalSpectrumGenerator` |
| Signal processing | `PeakPickerHiRes`, `GaussFilter`, `SavitzkyGolayFilter`, `Normalizer`, `ThresholdMower`, `WindowMower`, `NLargest`, `SpectrumAlignment`, `SpectrumAlignmentScore` |
| Features | `Feature`, `FeatureMap`, `ConsensusMap`, `MassTraceDetection`, `ElutionPeakDetection`, `FeatureFindingMetabo` |
| Identification | `PeptideHit` (with `getPeakAnnotations()`), `PeptideIdentification(List)`, `ProteinHit`, `ProteinIdentification(List)`, `SimpleSearchEngineAlgorithm`, `ProSEAlgorithm` (`search()`, `searchFiles()`), `ProSEResultSet`, `annotateSpectrum()`, `annotatePeaks()`, `predictMS2()` |
| Parameters | `Param`, and on every algorithm `getParameters()`/`setParameters()`, `getParameterValues()`/`setParameterValues({...})` |
| Misc | `getVersion()`, `getOpenMSDataPath()`, `setLogEnabled(channel, enabled)`, `getThreadingInfo()`, `getMemoryInfo()`, the file system `FS` with `mountHostDirectory()` and `unmount()`, files in parts with `allocateFile()`, `writeFileChunk()` and `readFileChunk()` |

`openms.d.ts` declares all of them for TypeScript, including the shapes of the plain objects
they return (e.g. `getPeaks(): { mz: Float64Array; intensity: Float32Array }`).

## Using the module

### Loading

`openms.mjs` is an ES module whose default export creates an instance of the
library. It loads `openms.wasm` from the directory it was loaded from; the file includes
OpenMS' shared data (modifications, enzymes, controlled vocabularies, ...), which the
library finds at `/share/OpenMS` of its file system. Serve `openms.wasm` so that browsers
can cache it (e.g. `Cache-Control: public, max-age=31536000, immutable` for a versioned
path): Chrome then also keeps the compiled code, and the module starts and runs faster from
the second visit on (see Performance below). Options of
[Emscripten's Module object](https://emscripten.org/docs/api_reference/module.html) can
be passed, e.g. to redirect OpenMS' log output or to locate the files elsewhere:

```js
const OpenMS = await createOpenMS({
  print: (line) => console.log(line),      // OpenMS' info log (stdout)
  printErr: (line) => console.warn(line),  // warnings and errors (stderr)
  locateFile: (file) => `/static/openms/${file}`,
});
OpenMS.setLogEnabled("info", false);       // silence progress messages
```

Browsers have to serve `openms.wasm` as `application/wasm`. Large files are best
processed in a Web Worker, so that the page stays responsive:

```js
// worker.mjs, started with new Worker("worker.mjs", { type: "module" })
import createOpenMS from "@openms/openms-wasm";
const ready = createOpenMS();                    // set onmessage before awaiting anything

onmessage = async ({ data: file }) => {          // a File from <input type="file">
  const OpenMS = await ready;
  OpenMS.FS.mkdir("/input");
  OpenMS.FS.writeFile(`/input/${file.name}`, new Uint8Array(await file.arrayBuffer()));
  const handler = new OpenMS.FileHandler();
  const exp = new OpenMS.MSExperiment();
  try {
    handler.loadExperiment(`/input/${file.name}`, exp);
    postMessage(exp.getSpectrumSummary());       // plain objects: no cleanup, cloneable
  } finally {
    handler.delete();
    exp.delete();
    OpenMS.FS.unlink(`/input/${file.name}`);
  }
};
```

`examples/node/summary.mjs` shows the same with Node.js, where the file's directory is
mounted instead of copied.

### Files

OpenMS reads and writes files by path. In WebAssembly, paths refer to the module's
[virtual file system](https://emscripten.org/docs/api_reference/Filesystem-API.html)
(WasmFS), `OpenMS.FS`, which keeps files in memory:

```js
// browser: a File from <input type="file"> or drag and drop
const bytes = new Uint8Array(await file.arrayBuffer());
OpenMS.FS.writeFile("/tmp/run.mzML", bytes);

const exp = new OpenMS.MSExperiment();
new OpenMS.MzMLFile().load("/tmp/run.mzML", exp);

// write a result and hand it to the user
new OpenMS.FeatureXMLFile().store("/tmp/features.featureXML", featureMap);
const result = OpenMS.FS.readFile("/tmp/features.featureXML");   // Uint8Array
```

`MzMLFile` also loads from and stores to strings (`loadBuffer(text, exp)`,
`storeBuffer(exp)`). In Node.js, `OpenMS.mountHostDirectory("/data", "/path/on/disk")`
makes a directory of the host visible at `/data`, without copying its files;
`OpenMS.unmount("/data")` removes the mount point again. Browsers have no such mount:
write the user's file into the file system (above). The file system is WasmFS, implemented
inside the module: every thread reads it directly, which matters for ProSE, which reads
the spectra on a helper thread while the other threads build the index. (With Emscripten's
JavaScript file system, which lives on the main thread, those reads waited for the busy
main thread; the index then could not be restricted to the spectra, and a search against
the human proteome took 8.7 s instead of 5.8 s.)

The module replaces Emscripten's open and close of host files (`wasmfs_node.js`; drop it when
[emscripten#27825](https://github.com/emscripten-core/emscripten/issues/27825) is fixed upstream):
Emscripten 6.0.12 took the error of a failed open, e.g. of a file deleted on the host after the
module had looked it up, for a file descriptor of the Node.js process, then read and closed it.
Errors are now reported as such. `FS.writeFile()` traps where its open of a host file fails (a
read-only file; it used to write to the process' stderr) and where its close fails (it used to go
on).

`FS.writeFile()` of WasmFS (Emscripten 6.0.12) appends to a file that exists, instead of replacing its
contents: `FS.unlink()` the file first, or use `allocateFile(path, size, capacity)`, which creates it
anew. (`allocateFile(path, size)` resizes an existing file in place, keeping its old bytes, so write all
of it.)

Gzip- and bzip2-compressed mzML files are read as usual.

**Files in parts.** Large files are written and read in chunks of a bounded size, so that JavaScript
never holds a whole file next to the module's copy:

```ts
allocateFile(path: string, size: number): void;                          // create (or resize) with size bytes, in one allocation
allocateFile(path: string, size: number, capacity: number): void;        // create anew with size bytes in a block of capacity bytes
writeFileChunk(path: string, data: Uint8Array, offset: number): void;    // write data at offset (the file is created if missing)
readFileChunk(path: string, offset: number, length: number): Uint8Array; // up to length bytes from offset (fewer at the end)
```

```js
// in a Web Worker: a File of the page, staged in chunks of 16 MB (FileReaderSync reads a slice synchronously)
OpenMS.allocateFile("/input/0/run.mzML", file.size);
const reader = new FileReaderSync();
for (let at = 0; at < file.size; at += 16 << 20) {
  OpenMS.writeFileChunk("/input/0/run.mzML", new Uint8Array(reader.readAsArrayBuffer(file.slice(at, at + (16 << 20)))), at);
}
```

Each chunk is copied once (450-850 MB/s in Chrome, with the reads of `FileReaderSync`; 7.7 GB/s writing
and 3.1 GB/s reading in Node.js). Allocate the file first when its size is known: the memory file
system keeps a file in one block, which grows by doubling otherwise and then briefly takes up to
three times its size. Errors, after which the module stays usable:

- A file has to be smaller than 2 GiB (at most 2080374784 bytes, 2 GiB - 64 MiB): a larger size or end
  of a write throws `InvalidValue: the value '<bytes>' was used but is not valid; a file of the module's
  file system has to be smaller than 2 GiB (at most 2080374784 bytes): <path>` before anything is
  allocated (an app can tell the user that the file is too large). A file of more than about 1 GiB
  cannot grow by a write past its end (that takes a block of 2 GiB): `InvalidValue: ... growing the file
  <path> of <n> bytes takes a block of ... allocate the file at its final size first (allocateFile())`.
- Where the module's memory cannot hold the file: `OutOfMemory`.
- `FS.writeFile()`, `FS.write()` and `FS.read()` throw an `FS.ErrnoError` of errno 48 (ENOMEM) there:
  Emscripten's took the 0 that `malloc` returns then for the address of their buffer and wrote over the
  start of the module's memory (`wasmfs_fs.js`, which also makes `FS.write()` copy a typed array of
  bytes in one piece instead of byte by byte). `FS.write()` and `FS.read()` throw errno 28 (EINVAL) for
  an `offset` and `length` beyond their buffer.

**Files one after the other: allocate them with a capacity.** The allocator (mimalloc) does not give the
block of a removed file to a larger one. Files that come one at a time, each removed before the next,
therefore grow the heap with each file that is larger than the ones before (0.9, 1.0 and 1.1 GiB took it
to 1103, 2147 and 3295 MB, and 1.2 GiB did not fit), unless their blocks have the same size.
`allocateFile(path, size, capacity)` gives the file a block of `capacity` bytes (at least `size`, smaller
than 2 GiB): allocate each file with the size of the largest as capacity, and every file takes the
block of the one before (0.9 to 1.4 GiB in increasing order: 1615 MB throughout). Only files allocated
this way, with the same capacity, reuse the block; the path may be the same or a new one. With a
capacity the file is created anew: an existing file at `path` is replaced (its contents go).

**Write within the size.** A write that goes past the size of a file makes WasmFS allocate a larger block
(twice the capacity) inside the file system, which is not guarded: where the memory cannot hold that
block, **a write past a file's capacity can abort the module** (a `RuntimeError`, fatal: see "What is not
guarded"). `writeFileChunk()` checks a write past the size as growth by doubling the size, but WasmFS does
not tell the capacity of a file, so it cannot know the block that a file with a capacity grows to. The
recipes allocate each file at its final size and write within it; do the same.

**What is not guarded.** These paths still allocate inside WasmFS without a check, and where the memory
cannot hold the file they abort the module: a write into the memory file system by `FS.writeFile()` or
`FS.write()` that grows the file (their own buffer is checked, the growth of the file is not),
`FS.truncate()`, and files that OpenMS writes, e.g. `storeIdXML()`, whose file grows by doubling (up to
three times its size) and fails when the heap is high. An exception that is not an `Error`, a raw
`WebAssembly.Exception` (e.g. `std::bad_alloc` from inside WasmFS) or a `RuntimeError` ("unreachable",
an abort), is fatal: the module's state is undefined (a lock of the file system may be held). Terminate
the Worker and start a new module; do not call the module again.

### Memory

Objects created from JavaScript (`new OpenMS.MSExperiment()`, but also objects returned
by functions, like `exp.getSpectrum(0)`) live in the WebAssembly memory and are not
garbage collected: call `.delete()` when done, or, where the JavaScript engine supports
[`using`](https://github.com/tc39/proposal-explicit-resource-management) declarations
(Node.js 24, current Chrome and Firefox), declare them with `using exp = new OpenMS.MSExperiment();`.
Plain values are converted to JavaScript and
need no cleanup: numbers, strings, arrays, peak data (`Float64Array`/`Float32Array`), and
the plain objects returned for small records (precursors, residues, modifications,
parameter entries, `toArray()`/`toObject()` views of results).

Getters return copies (as in pyOpenMS): to change a spectrum of an experiment,
`getSpectrum(i)`, modify it, and `setSpectrum(i, spectrum)`.

### Parameters

Algorithms take their settings from a `Param` object, as in C++. From JavaScript, the
simplest is to read and set them as plain objects:

```js
const picker = new OpenMS.PeakPickerHiRes();
picker.setParameterValues({ signal_to_noise: 2, "SignalToNoise:win_len": 300 });
const centroided = new OpenMS.MSExperiment();
picker.pickExperiment(profileExperiment, centroided);
```

`getDefaults().getEntries()` lists every parameter with its type, default, description,
allowed values and limits, e.g. to build a settings form. Invalid values are rejected
when the parameters are set.

### Database search with ProSE

`ProSEAlgorithm.search(spectra, database, proteins, peptides)` searches an mzML file against a
FASTA file as the ProSE tool does for one input file: the search with PSM-level FDR
(`FDR:PSM`), then, if `FDR:protein` is set and the database has decoys, protein inference with
picked-protein FDR. By default (parameter `rescore`) the search rescores the PSMs with Percolator
before its FDR, in-process as the tool does (needs decoys and at least 100 PSMs, else it is
skipped; the FDR filters then use its q-values); `search(..., { rescore: false })` sets the
parameter for the call and keeps ProSE's scores. It fills the two lists and returns a summary: the statistics of the run
and of the fragment index, Percolator's outcome (`applied`, why it was `skipped`, and a `warning` if its
sanity check replaced the learned model by the best single feature), and for open searches the
modification analysis (delta-mass histogram, mapped PTMs).

```js
const prose = new OpenMS.ProSEAlgorithm();
prose.setParameterValues({ "fragment:mass_tolerance": 0.02, "modifications:fixed": ["Carbamidomethyl (C)"] });
prose.setProgressCallback((p) => postMessage(p));   // { event, label, depth, begin, end, value, file }
const proteins = new OpenMS.ProteinIdentificationList(), peptides = new OpenMS.PeptideIdentificationList();
const summary = prose.search("/input/run.mzML", "/input/db.fasta", proteins, peptides);
summary.run.targetPSMs; summary.search.indexedPeptides; summary.modificationAnalysis?.ptms.entries;
```

The progress callback is called while the search runs, on the thread that set it: in a Web
Worker, forward the reports to the page with `postMessage()`. Each hit carries the fragment ions
ProSE matched (`hit.getPeakAnnotations()`, on the preprocessed spectrum); `annotateSpectrum(spectrum,
peptide, charge, { tolerance, unit, ions, losses })` matches the ions of a peptide to the peaks of
a spectrum as measured, for spectrum viewers. `getMemoryInfo().heapSize` after a search is its
memory high-water mark (WebAssembly memory does not shrink).

The ProSE app (<https://timosachsenberg.github.io/prose-web/>; source:
<https://github.com/timosachsenberg/prose-web>, the application in `app/`) is a web application
around it: see "Several files and the ProSE app" below. `examples/node/prose.mjs` runs a search
from the command line and reports time and memory.

Memory is the limit of a browser search: the module can address 4 GB, and a search that needs
more stops the module (a WebAssembly trap; start a new one). Measured with Node.js (threaded
build, 16 threads) on a 60-minute Orbitrap TMT run (6,103 MS2 spectra, 78 MB mzML):

| Database, search | Indexed peptides | Memory | Time | Native ProSE (memory, time) |
|---|---|---|---|---|
| *Erwinia* (4,499 proteins), tryptic | 0.5 M | 0.9 GB | 2.9 s | 0.25 GB, 2.5 s |
| *Erwinia*, open search (±500 Da) | 0.5 M | 0.9 GB | 3.9 s | 0.46 GB, 5.7 s |
| human (20,652 proteins), tryptic | 4.0 M | 2.0 GB | 8.6 s | 1.6 GB, 6.2 s |
| human, semi-tryptic | | more than 4 GB: stops | | 28 GB, 47 s |
| human, semi-tryptic, `database:chunk_size` 2000 | 77 M | 3.8 GB | 24 s | 1.2 GB, 43 s |

(Native: the ProSE tool of June 2026 in a container, which adds about a second of startup.)
`database:chunk_size` builds and searches the index in parts and so keeps larger searches under
4 GB; `low_memory` does the same for non-specific searches. A 64-bit module (Memory64; the
experimental presets `wasm64-emscripten-mt-release` and `wasm64-emscripten-mt-onnx-release`, see
"Building") raises the limit to 16 GB (Chrome, Firefox and Node.js; Safari does not support
Memory64). Measured the same way (October 2026), it runs the semi-tryptic human search without
chunking (54 s, 10 GB) and with `database:chunk_size` 8000 (36 s, 10.5 GB), which stop at 4 GB in
wasm32, with the same results; but it is 10 to 20% slower, and with Emscripten's allocator (mimalloc)
the module's memory grows far more there: 1.9 instead of 0.9 GB for *Erwinia*, 2.5 instead of 1.3 GB
for human, 12 instead of 2.5 GB for the semi-tryptic human search with `database:chunk_size` 2000
(resident: 1.1 and 1.3 GB, about as in wasm32, but 10.5 instead of 1.8 GB). More threads than about
8 to 16 hardly shortened these searches: decoy generation, parts of the index build and reading
the spectra are serial.

#### Several files and the ProSE app

The ProSE app (<https://timosachsenberg.github.io/prose-web/>, a web application for experiments
of several runs) is built on `searchFiles()` and the functions below.
`searchFiles(files, database, options)` searches several mzML files against one FASTA file and
pools the results as the ProSE tool's `-out_merged` does: the PSMs of all files, protein
inference across the files and, with `FDR:protein` and decoys in the database, the picked-protein
FDR. It returns a `ProSEResultSet`, which hands the results over as typed arrays:

```ts
searchFiles(files: string[], database: string, options?: ProSEFilesOptions): ProSEResultSet;
type ProSEFilesOptions = {
  rescore?: boolean; keepSpectra?: boolean; removeDecoys?: boolean;
  provideFile?: ((i: number, path: string) => void) | null;   // file by file (below)
  releaseFile?: ((i: number, path: string) => void) | null;
} | null | undefined;
interface ProSEResultSet {
  summary(): ProSEMultiSummary;       // per file: counts and Percolator's outcome; pooled: FDR levels, flags, modification analysis
  tables(): ResultTables;             // PSMs, peptides and proteins as columns
  spectra(): StoredSpectra;           // the peaks of the identified spectra
  storeIdXML(path: string): void;     // the pooled results as idXML
  delete(): void;
}
type ProSEMultiSummary = {
  files: ProSEFileSummary[];
  search: ProSESummary['search'];     // as search() reports it: database, enzyme, decoyMode, indexedPeptides, indexedFragments, secondsIndexBuild, secondsTotal, ...
  proteinFdrApplied: boolean; psmFdr: number; proteinFdr: number;
  decoyString: string; openSearch: boolean; modificationAnalysis: ProSESummary['modificationAnalysis'];
};
type ProSEFileSummary = {
  path: string; status: 'EXECUTION_OK' | 'INPUT_FILE_EMPTY' | 'UNEXPECTED_RESULT' | 'ILLEGAL_PARAMETERS' | 'UNKNOWN_ERROR';
  ms2Spectra: number; targetPSMs: number; decoyPSMs: number;
  percolator: { applied: boolean; skipped: string | null; warning: string | null } | null;
  peptdeepNce: number | null;         // the collision energy the PeptDeep predictions of the file used
  precursorError: ProSESummary['run']['precursorError'];   // { median, mad, recommended } of the PSMs of the file; null: not estimated
  fragmentError: ProSESummary['run']['fragmentError'];
  secondsSearch: number;
};
```

- `rescore` sets the parameter `rescore` for the call, as for `search()` (default: the parameter, true).
- `keepSpectra` (default true) keeps the peaks of the identified spectra for `spectra()`.
- `removeDecoys` (default false): with the protein FDR, remove the decoy PSMs and proteins from the
  pooled results, as `-out_merged` does. By default the decoys that pass are kept, to show them.
- The PSMs name their file by its path, so repeated paths are rejected. Files of the same name go
  to different directories (`/input/0/run.mzML`, `/input/1/run.mzML`).
- If no file has a PSM, the result is empty (no error): no PSMs, peptides or proteins (`targetPSMs` and
  `decoyPSMs` are 0 and the tables have count 0; `ms2Spectra` still counts the file's spectra), no
  protein FDR, an empty idXML.
- Percolator needs decoys and 100 PSMs per file. If it rescored only some of the files, their
  scores cannot be pooled with the others', and `searchFiles()` throws (`InvalidParameter`, naming the
  files that were not rescored): search with `{ rescore: false }`, or leave those files out.

- `provideFile` and `releaseFile`: search the files one at a time, each while it is there (below).
- `summary().search` holds the statistics `search()` has under `search`, for the whole search: the database,
  the enzyme, the decoy mode (`decoyMode`: `generated`, `external` or `none (target-only)`), the size of the
  fragment index (`indexedPeptides`, `indexedFragments`) and the times. The index is built once for all
  files, so `secondsIndexBuild` is its time; `secondsTotal` is the time of the whole `searchFiles()` call,
  with pooling and reading the spectra. (The index of several files holds all peptides of the digest. That
  of `search()`, and of a search of one file, holds only those its spectra can reach when the search runs
  on one thread; with more threads it, too, holds the whole digest. So the sizes differ with the number of
  files and of threads.) With
  `provideFile` the index is built for each file, from its spectra: `secondsIndexBuild` adds up the builds,
  the size is that of the first file's index, and `secondsTotal` includes the writing and removing of the
  files. The tolerances ProSE recommends from the mass errors of the PSMs (`precursorError`,
  `fragmentError`: `{ median, mad, recommended }` or `null`) are those of each file, so they are in
  `files[i]`, not in `search`.

**Inputs larger than the memory: file by file.** The input files are searched while they are in the
module's file system, which in a browser keeps them in the module's memory (4 GB at most, and
WebAssembly memory never shrinks): with all files written first, two files of 733 MB already take
the heap to 3.9 GB. With `provideFile`, `searchFiles()` searches the files one at a time:

1. `provideFile(i, path)` writes file `i` (its position in `files`) at `path`, e.g. from the user's
   `File` with `allocateFile()` and `writeFileChunk()` (above). It is called synchronously from the
   search and has to have written the file when it returns, so in a browser it reads the `File` with
   `FileReaderSync`: run the module in a Web Worker (as recommended anyway). It must be synchronous:
   an `async` function, or one that returns a Promise (any thenable), makes `searchFiles()` throw
   `IllegalArgument: searchFiles(): provideFile must be synchronous (return nothing, not a Promise):
   provideFile(0, '<path>') returned a Promise.` (the same for `releaseFile`).
2. The file is searched on its own, and the spectra of its PSMs are read for `spectra()`.
3. `releaseFile(i, path)` removes it, e.g. with `OpenMS.FS.unlink(path)`.

```js
// in a Web Worker: files are the user's File objects (from the page), the database is written already
const paths = files.map((f, i) => `/input/${i}/${f.name}`);   // distinct paths, as always
const chunk = 16 << 20, reader = new FileReaderSync();
// one block of the size of the largest file for all of them (allocateFile with a capacity, above): files of
// increasing size would grow the heap with every file otherwise. Larger than 2080374784 bytes: too large.
const capacity = Math.max(...files.map((f) => f.size));
const set = prose.searchFiles(paths, "/input/database.fasta", {
  provideFile: (i, path) => {
    OpenMS.FS.mkdirTree(path.slice(0, path.lastIndexOf("/")));
    OpenMS.allocateFile(path, files[i].size, capacity);
    for (let at = 0; at < files[i].size; at += chunk) {
      OpenMS.writeFileChunk(path, new Uint8Array(reader.readAsArrayBuffer(files[i].slice(at, at + chunk))), at);
    }
  },
  releaseFile: (i, path) => { try { OpenMS.FS.unlink(path); } catch { /* not written */ } },
});
```

Then the results of all files are pooled as without `provideFile`. Only the largest file has to fit
into the memory, next to the results. The results are the same as with all files there (the tests
compare them, also for open searches). What differs:

- The fragment index is built for each file (from its spectra), instead of once for all files, and
  the next file is not read while one is searched: in Chrome with 4 threads, 13.9 s instead of 14.5 s
  for one file of 733 MB (plus 1.2 s to write it), and about 13 s for each further file.
- Progress: each file has its step "Searching <name> (i of n)..." as before, with `file` = i - 1 for it
  and everything in it, writing and removing the file included; "Reading the spectra of the PSMs..."
  comes in the step of each file (depth 1, with its `file`), not once after the files.
- What `provideFile` or `releaseFile` throws ends the search: `searchFiles()` throws an `Error` naming
  the function, the file and its message (`SearchFilesError: searchFiles(): provideFile(1, '/input/1/run.mzML')
  threw: ...`). `releaseFile` is also called for the file whose writing or search failed (what it throws
  then is ignored), so a `releaseFile` that removes the file if it is there cleans up in every case.
- `ProSEResultSet._debugSetNativeId()` reads the input files again, which are gone.

Without `provideFile`, all files have to be there for the whole call: in Node.js, mount their directory
(`mountHostDirectory()`, no copy); in a browser, write them first, which is fine for a total well below
the memory (a few hundred MB).

Measured in Chromium 141 (Playwright), 4 threads, a closed search with Percolator of copies of a 733 MB
mzML file (48824 MS2 spectra, 10011 target PSMs each); heap is `getMemoryInfo().heapSize` after the search,
the time that of `searchFiles()` without writing the files:

| Inputs | All files written first | File by file |
|---|---|---|
| 1 × 733 MB | 2393 MB², 14.5 s | 2153 MB, 13.9 s (+ 1.2 s writing) |
| 2 × 733 MB | 3969 MB¹, 23.6 s | 2269 MB, 27.0 s (+ 2.4 s) |
| 4 × 733 MB | does not fit | 2485 MB, 52.9 s (+ 5.1 s) |
| 12 × 733 MB | does not fit | 2969 MB, 163 s (+ 12 s) |

¹ ² written with `FS.write()` in chunks of 16 MB, by which a file grows by doubling (`allocateFile()` first
avoids that). ¹ is from a run before this module's `FS.write()`, which copied byte by byte (the same run
gave 2569 MB and 13.7 s for one file); ² and the file-by-file column are from the module as it is.

The heap still grows with the number of files, by about 70 MB per file of this size (110 MB for the
second): part of it are the results (the PSMs and the peaks of their spectra), the rest is memory the
allocator (mimalloc, with threads) does not reuse. With many small files it is more, and it grows with
the number of threads and with PeptDeep: 12 files of 75 MB, file by file, ended at 3007 MB with 8 threads
and at 3971 MB with 16 (PeptDeep off), and at 3263 MB with 4 threads, 4079 MB with 8 and the 4 GB limit
with 16 (PeptDeep on). For searches of many files, use fewer threads (4 kept 12 such files with PeptDeep
at 3.3 GB). A file has to be smaller than 2 GiB (above).

**Exporting idXML.** `storeIdXML(path)` writes into the memory file system as well, where the file
grows by doubling and briefly takes up to three times its size (not guarded: see "Files in parts"). Hand
it over in chunks and remove it:

```js
set.storeIdXML("/output/merged.idXML");
const size = OpenMS.FS.stat("/output/merged.idXML").size, parts = [];
for (let at = 0; at < size; at += 16 << 20) parts.push(OpenMS.readFileChunk("/output/merged.idXML", at, 16 << 20));
OpenMS.FS.unlink("/output/merged.idXML");
const blob = new Blob(parts, { type: "application/xml" });
```

Call `delete()` on the result set when done: its results, among them the heap copy of the peaks,
stay in the module's memory until then, and WebAssembly memory never shrinks.
(`_debugSetNativeId(row, nativeId)` is for the tests only: it replaces the native ID of a PSM and
reads the spectra from the input files again.)

**Progress.** The reports of `setProgressCallback()` have one more field, `file`. Each input file is
a step of its own, "Searching <name> (i of n)...", with `file` = i - 1 for it and the steps nested
in it (<name> is the file name without its directory). After the files comes "Reading the spectra
of the PSMs..." with `file` = -1 (with `keepSpectra`, if there are PSMs). A chunked search (`database:chunk_size`) has no
per-file step: `file` is -1 throughout.

**`tables()`** returns the pooled results as typed columns (copies, which JavaScript owns):

```ts
type ResultTables = {
  files: string[];                    // the input files, in the order of the search
  psms: {                             // one row per identified spectrum (its best hit)
    count: number;
    file: Uint16Array;                // position in files
    nativeId: string[]; rt: Float64Array; mz: Float64Array; charge: Int8Array;
    peptide: Uint32Array;             // row of peptides
    score: Float64Array; qValue: Float64Array; pep: Float64Array;
    precursorErrorPpm: Float64Array; isotopeError: Int8Array;
    deltaMass: Float64Array;          // mass shift of an open search; NaN: none
    decoy: Uint8Array;
    featureNames: string[];
    features: Float64Array;           // row-major, count * featureNames.length; NaN: absent
  };
  peptides: {                         // unique by openms
    count: number;
    openms: string[]; proforma: string[]; sequence: string[]; decoy: Uint8Array;
    modOffsets: Uint32Array;          // modifications of peptide k: modOffsets[k] .. modOffsets[k + 1] of
    modPosition: Int16Array; modName: string[]; modUnimod: Int32Array; modMass: Float64Array;   // these
    proteinOffsets: Uint32Array;      // proteins of peptide k: proteinOffsets[k] .. proteinOffsets[k + 1] of
    protein: Uint32Array; start: Uint32Array; end: Uint32Array;                                   // these
  };
  proteins: {
    count: number;
    accession: string[]; description: string[]; sequence: string[];
    score: Float64Array; qValue: Float64Array; decoy: Uint8Array;
    group: Int32Array;                // index of the indistinguishable-protein group; -1: none
  };
};
```

- `psms.score` is Percolator's (SVM) score if the search rescored the PSMs, else ln(hyperscore); higher
  is better. `proteins.score` is the best `psms.score` of the protein's PSMs. `pep` is NaN where the file
  was not rescored. `psms.qValue` is Percolator's q-value for a rescored file; otherwise it is the q-value of
  the target-decoy competition of the file's PSMs, as `FDR:PSM` computes it, **also with `FDR:PSM` 0**, which
  filters nothing (the scores stay ProSE's: only the q-values are added). It is NaN only when the search had no
  decoys (or the file no decoy PSM). `proteins.qValue` is NaN without a protein FDR.
- `featureNames` holds the features Percolator trained on, those that some PSM has: `ln(hyperscore)`,
  the standard features of the .pin format (`CalcMass`, `mass`, `peplen`, `charge<N>`, `enzN`, `enzC`,
  `enzInt`, `dm`, `absdm`, computed also for searches without rescoring; the enzyme is
  `ProSEAlgorithm::pinEnzymeName()` of the search's), the search engine's extra features and PeptDeep's.
- A peptide `k` with `proteinOffsets[k] === proteinOffsets[k + 1]` has no protein: all its proteins failed the
  protein FDR, and the references to them were removed. The peptide and its PSMs stay (19 of 967
  peptides in one test search), so an application can show them as "below protein FDR" and leave them out
  of the protein counts. Without a protein FDR every peptide has proteins.
- `peptides.proforma` is ProForma with Unimod names (`AVM[Oxidation]DDFAAFVEK`). It never has the
  unlocalized mass shift of an open search: that is a value of each PSM. Compose it from
  `psms.deltaMass`, which is NaN in a closed search. In an open search the shift of an unmodified
  PSM is its precursor mass error (typically a few mDa); ProSE counts |deltaMass| <= 0.05 Da as
  unmodified (`unmodifiedPSMs` of `summary().modificationAnalysis.deltaMasses`), so test for a shift
  above that:

  ```js
  const d = psms.deltaMass[i];
  const text = Number.isFinite(d) && Math.abs(d) > 0.05
    ? `[${d < 0 ? "-" : "+"}${Math.abs(d).toFixed(4)}]?${proforma}` : proforma;   // e.g. [-57.0215]?PEPTIDEK
  ```
- `proteins.sequence` is the sequence as searched: from the FASTA, for a decoy of the FASTA too. A decoy
  that ProSE generated has the sequence it generated (regenerated with ProSE's rule), so that
  `peptides.start`/`end` (0-based, end exclusive) slice the peptide from it. `''` where unknown.

**`spectra()`** returns the MS2 peaks of the identified spectra, read from the input files at the end
of the search (`keepSpectra`):

```ts
type StoredSpectra = {
  psmSpectrum: Int32Array;            // per row of tables().psms: the spectrum; -1: none
  offsets: Uint32Array;               // spectrum s: offsets[s] .. offsets[s + 1] of mz and intensity
  mz: Float32Array; intensity: Float32Array;
  precursorMz: Float64Array;          // per spectrum; NaN: none
  collisionEnergy: Float32Array;      // per spectrum, in eV; NaN: unknown
};
```

`collisionEnergy` is the collision energy MS:1000045 that converters write (msconvert, ThermoRawFileParser, Bruker),
which OpenMS keeps as the precursor's meta value "collision energy" (PeptDeep's rescoring reads that one); where a
file has none, the activation energy MS:1000509. NaN where it has neither (OpenMS has 0 for an unknown activation energy).

The peaks are the raw MS2 peaks of the file (not preprocessed), as Float32 arrays sorted by m/z.
A row keeps its spectrum if `qValue <= psmFdr` (decoys included), if `FDR:PSM` is 0 (no filter:
every row), or if `qValue` is NaN. A PSM has -1 where `keepSpectra` is false, where its `qValue`
is above `psmFdr`, where the file has no MS2 spectrum with its native ID, or where the file could
not be read again. Where a file has
several spectra with the same native ID, the PSM gets the one at its retention time (within a
microsecond), and -1 if not exactly one is. PSMs of one spectrum share its index.

**`annotatePeaks(mz, intensity, peptide, charge, options)`** matches the fragment ions of a peptide
to peaks given as arrays, e.g. `subarray()`s of the stored ones (`annotateSpectrum()` takes an
`MSSpectrum`; the options and the matching are the same). `options` is required by the declaration
(the module, built without assertions, does not check the number of arguments): pass `null` for the
defaults, a tolerance of 0.02 `Da`, `maxCharge` of the precursor charge - 1 (at least 1), `ions`
`['b', 'y']`, no `losses` and no `precursor` peaks:

```ts
annotatePeaks(mz: ArrayLike<number>, intensity: ArrayLike<number>, peptide: string, charge: number, options: AnnotationOptions): PeakAnnotationResult;
type AnnotationOptions = {
  tolerance?: number; unit?: 'Da' | 'ppm'; maxCharge?: number;
  ions?: ('a' | 'b' | 'c' | 'x' | 'y' | 'z')[]; losses?: boolean; precursor?: boolean;
} | null | undefined;
type PeakAnnotationResult = {
  matches: PeakMatch[];               // by ascending m/z
  theoretical: { ion: string; type: string; ordinal: number; charge: number; mz: number }[];   // every ion, matched or not
};
type PeakMatch = {
  peak: number;                       // index in mz and intensity as passed
  ion: string; type: string; ordinal: number; charge: number;
  loss: string;                       // formula of the neutral loss; '' if none
  theoreticalMz: number;
  errorPpm: number;                   // (observed - theoretical) / theoretical * 1e6, whatever the unit
};
```

- In `ppm`, every ion takes its closest peak within the tolerance, so a peak can carry several ions
  (ions of one m/z, e.g. `z1++` and `y1-NH3++`): several matches have the same `peak`. Of peaks of
  equal m/z, the first in the arrays is taken.
- In `Da`, the matching is one to one in m/z order: a peak has at most one ion (of two ions of one m/z,
  one is left out, unless there are two peaks). Of peaks of equal m/z one is taken; which one is
  unspecified.
- `maxCharge` defaults to the precursor charge - 1 (at least 1): a precursor of charge 3 gets b and y ions of
  charge 1 and 2. `predictMS2()` gives doubly charged ions for a precursor of charge **2** or more, so
  the two disagree for charge 2: `annotatePeaks()` then has singly charged ions only. Pass `maxCharge`
  explicitly where they have to agree (`maxCharge: 2` for the mirror plot of `predictMS2()`; the ProSE app does).
- Ions are named `b3`, `y5++` (one `+` per charge, none for charge 1), `y7-H2O++` with a neutral loss,
  `[M+2H]++` and `[M+H-NH3]` for precursor peaks. The losses are written as formulas: `H2O`, `NH3`,
  `H3PO4` (phosphorylation) and `HPO3`; any other as the empirical formula without counts of 1
  (`CH4OS` of oxidized methionine). Precursor peaks are generated for the charge `maxCharge` only.

**`predictMS2(peptide, charge, nce, instrument, modelPath?)`** predicts the fragment intensities of a
peptide with the PeptDeep MS2 model, for a spectrum viewer's mirror plot (ONNX builds only; a module
without ONNX Runtime throws "predictMS2: this OpenMS module was built without ONNX Runtime"):

```ts
predictMS2(peptide: string, charge: number, nce: number, instrument: string, modelPath?: string | null): PredictedIon[];
type PredictedIon = {
  ion: string;                        // as annotatePeaks names it: 'b3', 'y5++'
  type: 'b' | 'y'; ordinal: number;
  charge: 1 | 2;                      // 2 only if the precursor has charge 2 or more
  mz: number;                         // from the theoretical spectrum, with the modifications as given
  intensity: number;
};
```

```js
OpenMS.predictMS2("AVM(Oxidation)DDFAAFVEK", 2, 30, "QE");   // 13 ions by ascending m/z: [{ ion: "y1", type: "y", ordinal: 1, charge: 1, mz: 147.1128, intensity: 0.628 }, { ion: "b2", ... }, ...]
```

- `nce` is the normalized collision energy as `peptdeep:nce` takes it (e.g. 30; the file summary's
  `peptdeepNce` is the one a search used). The instrument is one of `QE`, `Lumos`, `timsTOF`,
  `SciexTOF`, `ThermoTOF` (case-insensitive); any other name goes to peptdeep's catch-all slot.
- The result is sorted by ascending m/z. Intensities are relative to the most intense returned b/y ion
  (which is 1), and ions of 0.01 or less are left out. A peptide of one residue gives `[]`.
- Ions are of charge 1, and also of charge 2 for a precursor of charge 2 or more (`annotatePeaks()` has
  charge 2 only from a precursor of charge 3, unless `maxCharge` is passed: see above).
- The charge has to be 1 to 6 and `nce` 10 to 60, else `InvalidValue` (`The precursor charge has to be 1 to 6 ...`,
  `The normalised collision energy has to be a number from 10 to 60 ...`). This range is the binding's choice, as
  neither the library nor peptdeep states one: the model answers for any value, but badly outside the usual ones
  (all charges from 8 on give the same ions, and an energy of 100 gives none). A search can use an energy
  outside it (`peptdeep:nce` takes any value, and the automatic one follows the energy in the file): clamp or
  show a message for such a file.
- `modelPath` defaults to `/models/peptdeep_ms2_dynamic.onnx` (write the model there, as described
  below). The ONNX session is cached per model path and never evicted (tens of MB per path; a path
  that later names another file keeps answering from the first one). A failed load is not cached.
- The first call takes about 0.5 s here (it starts ONNX Runtime's global thread pool and loads the
  model), the next few tens of ms (up to about 0.1 s for a new peptide length), later calls about
  3 ms: call it a few times early to hide the start.

### Machine learning: ONNX Runtime and PeptDeep

A module built with the CMake preset `wasm32-emscripten-mt-onnx-release` (`WITH_ONNX`, threaded)
contains ONNX Runtime; `OpenMS.hasONNXRuntime()` says whether a module has it. ProSE's
`peptdeep:enable` then works as natively: PeptDeep (AlphaPeptDeep) predicts the fragment
intensities and the retention time of every PSM, and their agreement with the spectrum becomes
five more Percolator features (`ms2_cosine`, `ms2_spectral_angle`, `ms2_pearson`,
`ms2_frac_pred_found`, `rt_abs_error`). The models are not part of the module: the build copies
them to `models/` of the package (MS2 15.5 MB, retention time 2.7 MB; from
https://archive.openms.de/openms/models). Write them to the module's file system and name them by
absolute path (`/share/OpenMS` of the module is read-only, so the default names do not resolve):

```js
const models = { "peptdeep:ms2_model": "peptdeep_ms2_dynamic.onnx", "peptdeep:rt_model": "peptdeep_rt_dynamic.onnx" };
OpenMS.FS.mkdir("/models");                         // Node.js: OpenMS.mountHostDirectory("/models", "<package>/models")
const parameters = { "peptdeep:enable": true };
for (const [parameter, name] of Object.entries(models)) {
  OpenMS.FS.writeFile(`/models/${name}`, new Uint8Array(await (await fetch(`models/${name}`)).arrayBuffer()));
  parameters[parameter] = `/models/${name}`;
}
prose.setParameterValues(parameters);
```

The ProSE app offers this as an option next to Percolator and loads the models on first use;
`examples/node/prose.mjs --peptdeep` mounts them. ONNX Runtime computes with WebAssembly SIMD on
one pool of as many threads as OpenMP, shared by all sessions (ONNX Runtime's WebAssembly build
uses global thread pools); `threads.js` starts their Web Workers with the module. The predictions are those of native OpenMS: the PeptDeep class tests,
compiled to WebAssembly, pass (largest difference from Python's ONNX Runtime: 1.2e-7 for retention
times, 1.1e-6 for fragment intensities), and with the same NCE the five features of 896 PSMs differ
from native ProSE's by at most 3.3e-7. Measured on the data of TOPP_ProSE_11_rescore (2,789 MS2
spectra, 896 PSMs, 16 threads):

| | Native ProSE | Node.js | Chromium |
|---|---|---|---|
| search with Percolator | 9.4 s | 10.6 s | 11.2 s |
| with PeptDeep features | 12.2 s | 13.9 s | 12.7 s (14.5 s the first time: models) |
| memory with PeptDeep | 0.17 GB | 0.8 GB | 1.0 GB |

ONNX Runtime makes `openms.wasm` 42.6 MB instead of 29.1 MB (gzip: 10.3 instead of 6.8 MB); the
models are 17 MB (gzip). The automatic NCE (`peptdeep:nce` -1) picks the grid point with the best
median cosine; where the grid is nearly flat, rounding differences between platforms can pick
another point (here 26 natively, 36 in WebAssembly, median cosines 0.55504 and 0.55503) and with it
other feature values.

### Errors

An OpenMS exception becomes a JavaScript `Error` whose message names it, e.g.
`ParseError: ... while parsing 'PEPT[X'`.

### Threads

The package holds the module twice, with the same API: `st/` single-threaded (the parts of
OpenMS parallelized with OpenMP run sequentially) and `mt/` with threads (OpenMP regions, e.g.
mzML decoding, feature finding, database search, and the helper threads of OpenMS run on
WebAssembly threads: Web Workers, in Node.js worker threads, that share the module's memory).
`createOpenMS()` starts the threaded build where threads are available:

- **Node.js**: always, nothing to set up.
- **Browsers**: on cross-origin isolated pages, i.e. served with the headers
  `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`
  (`credentialless` also works in Chrome and Firefox, not in Safari): only they may share
  memory between threads. Elsewhere the single-threaded build.

Which one you got, and with how many threads:

```js
import createOpenMS, { threadSupport } from "@openms/openms-wasm";

threadSupport();   // before loading: { supported: false, reason: "the page is not cross-origin isolated (...)" }
const OpenMS = await createOpenMS({ numThreads: 4 });   // numThreads is optional
OpenMS.getThreadingInfo();   // { threaded: false, available: false, numThreads: 1 }
```

- `getThreadingInfo()`: `threaded`, the build that runs; `available`, whether it can start
  threads; `numThreads`, the threads OpenMS computes with.
- `threadSupport()` gives the same decision before loading, with the reason.
- Asking for more than one thread (`numThreads`) where there are none logs one console
  warning with the reason.
- `createOpenMS({ threads: false })` (or `true`) forces a build; `@openms/openms-wasm/st` and
  `@openms/openms-wasm/mt` import one directly, without the loader. (The threaded build also
  runs without shared memory, on one thread.)

`numThreads` defaults to `OMP_NUM_THREADS` of the process, else all logical cores, in
Node.js, and to half the logical cores, at most 4, in browsers, where the Web Workers for
these threads are started with the module. Load the module in a Web Worker: on the page's
main thread, waiting threads busy-wait and block the page. In Chromium, the database search
below took 1832 ms on 1 thread, 1342 ms on 2 and 1233 ms on 4 (in a Worker, on an isolated
page); Firefox and Safari are not tested yet.

What threads gain depends on how much of a workflow OpenMS parallelizes. Measured with
Node.js on 4 cores (best of 3 runs): a database search (SimpleSearchEngine, 1684 spectra of
a BSA run against an E. coli target-decoy database) took 1458 ms on 1 thread, 1013 ms on 2
and 880 ms on 4 (1549 ms in the single-threaded build); loading that mzML file (about 200 ms)
and the feature finding of FeatureFinderMetabo_3 (about 220 ms) took the same time on 1 to
4 threads.

### Serving the module

- `openms.wasm` has to be served as `application/wasm` (for `WebAssembly.instantiateStreaming`).
- Threads need a cross-origin isolated page: the headers `Cross-Origin-Opener-Policy: same-origin`
  and `Cross-Origin-Embedder-Policy: require-corp` on the page and on the scripts of the Web
  Workers it starts. (The module's own threads need none: in browsers their Workers are started
  from a Blob that imports the module, and so inherit the isolation of their creator.)
- Let browsers cache `openms.wasm` and `openms.mjs` (e.g. `Cache-Control: public, max-age=31536000,
  immutable` for versioned paths). Chrome keeps the compiled code of a module it could cache
  and skips the compilation on later visits; the module then also runs faster: V8 compiles
  WebAssembly in tiers, and code that was optimized once comes back optimized. A search against
  the human proteome (16 threads) took 5.9 s on the first visit and 3.7 s on the next; natively
  2.4 s. With `Cache-Control: no-store`, every visit is the first.
- Where the server does neither, a service worker can: the demo application (`examples/browser/`)
  registers `examples/sw.js` (through `examples/sw-client.mjs`), which adds the isolation headers
  to the page and to worker scripts, so threads work on static hosts such as GitHub Pages, and stores
  the module files in the Cache API, answering even the first request from the store, so that
  Chrome keeps the compiled code from the first visit on. On the first visit the page reloads
  once. The stores are named after the package version; a hard reload (Shift+Reload) bypasses
  the worker, e.g. after rebuilding the module with the same version.

### Limits

The module is a 32-bit (wasm32) program: it can address up to 4 GB (the experimental 64-bit
builds up to 16 GB, see "Building"). It cannot start
other programs (search engine adapters, R, Java), has no network access of its own
(`NetworkGetRequest` reports an error; fetch data in JavaScript and write it to the file
system), and cannot read vendor formats (Thermo RAW, Bruker .d): convert them to mzML
first. The HMDB mapping tables of accurate mass search are not packaged (write them to
`/share/OpenMS/CHEMISTRY/` if needed). A JavaScript engine reports no CPU time: the CPU
times of `StopWatch` are the elapsed times of the module.

## Building

The module is built by the OpenMS build system with the
[Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html) (tested with
Emscripten 6.0.12). vcpkg builds the dependencies for the `wasm32-emscripten` triplet
(`vcpkg-overlays/triplets/wasm32-emscripten.cmake`), which takes about an hour the first
time; building them needs `flex` and `bison` on the build machine (Apache Thrift).

```bash
git submodule update --init vcpkg
source /path/to/emsdk/emsdk_env.sh           # sets EMSDK, puts emcc on PATH
cmake --preset wasm32-emscripten-release      # build/wasm32-emscripten-release
cmake --build build/wasm32-emscripten-release --target openms_wasm -j$(nproc)
ctest --test-dir build/wasm32-emscripten-release -R openms_wasm --output-on-failure
```

The threaded build comes from the preset `wasm32-emscripten-mt-release` (build directory
`build/wasm32-emscripten-mt-release`): it sets `OPENMS_WASM_THREADS=ON` and builds the
dependencies a second time, for the `wasm32-emscripten-mt` triplet (with `-pthread`), as a
module with threads can only be linked from objects compiled for it. Each build writes a
package of its own module (`@openms/openms-wasm`, `@openms/openms-wasm-mt`); the published
package, with both builds and the loader, is assembled from the two:

```bash
cmake -DSINGLE=build/wasm32-emscripten-release/openms-wasm \
      -DMULTI=build/wasm32-emscripten-mt-release/openms-wasm \
      -DOUTPUT=build/openms-wasm-package -P src/openms_wasm/assemble_package.cmake
```

The preset `wasm32-emscripten-mt-onnx-release` is the threaded build with ONNX Runtime (`WITH_ONNX`,
vcpkg feature `onnxruntime`; see "Machine learning" above). vcpkg builds ONNX Runtime for the
`wasm32-emscripten-mt` triplet in about two minutes, with the ports of
`vcpkg-overlays/emscripten/ports`: `onnxruntime` as static libraries with WebAssembly SIMD and
threads, without LTO, onnxruntime-web's JavaScript module or cpuinfo; `onnx` without static schema
registration (ONNX Runtime registers the schemas itself). The build downloads the PeptDeep models
and copies them to `models/` of the package.

The presets `wasm64-emscripten-mt-release` and `wasm64-emscripten-mt-onnx-release` are the threaded
builds, without and with ONNX Runtime, for 64-bit memory (Memory64): the module's memory grows up to
16 GB instead of 4 GB (Chrome, Firefox and Node.js; not Safari). vcpkg builds the dependencies, ONNX
Runtime included, for the `wasm64-emscripten-mt` triplet (`-m64`; its chainloaded toolchain
`vcpkg-overlays/toolchains/wasm64-emscripten.cmake` sets the pointer size). They are experimental:
Embind passes a `size_t` as BigInt there, so `AASequence.size()`, `MSExperiment.size()`,
`getNrSpectra()`, `getNrChromatograms()` and `CoarseIsotopePatternGenerator.getMaxIsotope()` return a
`bigint`, and `openms.d.ts` declares indices and counts such as the one of `getSpectrum()` as `bigint`
(numbers are accepted). Emscripten's allocator (mimalloc on top of emmalloc, Emscripten 6.0.12) cannot
allocate a block of 4 GB or more in such a module, and takes more memory than in wasm32 (see the
ProSE measurements above).

The npm package is written to `build/wasm32-emscripten-release/openms-wasm/`
(`openms.mjs`, `openms.wasm`, `openms.d.ts`, `package.json`); `npm pack`
there creates the tarball. `examples/browser/` in it is a demo web application: serve the
package directory (`python3 -m http.server -d build/wasm32-emscripten-release/openms-wasm`)
and open `http://localhost:8000/examples/browser/`.

`cmake/emscripten.cmake` holds the settings of the WebAssembly build: static libraries;
no GUI, TOPP tools, OpenMP, curl, Thermo RAW or Bruker support; C++ exceptions as
WebAssembly exceptions (`-fwasm-exceptions`; all code of a module has to use the same, so
the triplet compiles every dependency with it through its chainloaded toolchain
`vcpkg-overlays/toolchains/wasm32-emscripten.cmake`); WebAssembly SIMD (`OPENMS_WASM_SIMD`); an 8 MB stack
(`OPENMS_WASM_STACK_SIZE`). Options of `src/openms_wasm/CMakeLists.txt` choose the shared
data that is packaged (`OPENMS_WASM_SHARE_DIRECTORIES`, `OPENMS_WASM_SHARE_EXCLUDE`) and
how: inside `openms.wasm` (`OPENMS_WASM_SHARE=EMBED`, the default), as separate
`openms.data` downloaded at start-up (`PRELOAD`; it is looked up relative to the page, so
pass `locateFile`), or not at all (`NONE`).

The OpenMS class tests can be run on the WebAssembly library as well: configure with
`-DENABLE_CLASS_TESTING=ON`, build the tests you want (each links the whole library) and
run them with `ctest`; they run in Node.js and use the files of the build machine. All of
them pass except those that exercise what the module does not have: other programs and the
network (`ExternalProcess`, `TOPPExternalToolBase`, `Network`; `File` looks for programs next
to its own), the TOPP tool registry (`ToolHandler`: it lists the tools of the build, and the
WebAssembly build has none), and, in the single-threaded build, threads
(`GlobalExceptionHandler`, `LogConfigHandler`, `LogStream`, `UniqueIdGenerator`, `File`,
`TOPPBase`: they start `std::thread`s; they pass in the multi-threaded build).

`tests/browser/run.mjs` tests the module in Chromium (with Playwright, as CI does): on
pages with and without cross-origin isolation, on the page and in a Web Worker, it checks
that the module starts, the threads it reports, the results of feature finding and a database
search, and that a Worker keeps the page responsive; `--benchmark` adds a larger search at
1, 2 and 4 threads.

```bash
npm install --prefix /tmp/pw playwright-core && node /tmp/pw/node_modules/playwright-core/cli.js install chromium
PLAYWRIGHT_MODULE=/tmp/pw/node_modules/playwright-core node tests/browser/run.mjs build/wasm32-emscripten-mt-release/openms-wasm
```

A parallel region of the multi-threaded build may capture at most 32 variables (a limit
of LLVM's OpenMP runtime on WebAssembly); a larger one exits the module ("Too many args to
microtask"). Wrapping the body of such a loop into a lambda keeps it at one.

The WebAssembly build installs no CMake package (`OpenMSConfig.cmake`): its libraries are
static and privately link bundled libraries that are not installed. To compile C++ code of
your own against OpenMS, add it to the module: `-DOPENMS_WASM_EXTRA_SOURCES=/path/my_bindings.cpp`
compiles your sources (with their `EMSCRIPTEN_BINDINGS`) into `openms.wasm`.

### Adding bindings

The bindings are [Embind](https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html)
declarations in `bindings/bind_<area>.cpp`. Wrap every bound function in `OMS_FN(...)`
(`bindings/binding_utils.h`): it turns C++ exceptions into JavaScript errors. Return
small records as plain objects (`val`) and peak data as typed arrays (`toTypedArray`);
bind classes the user creates, passes to algorithms or keeps. Add a test to `tests/`.
