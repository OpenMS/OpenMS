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
| Identification | `PeptideHit` (with `getPeakAnnotations()`), `PeptideIdentification(List)`, `ProteinHit`, `ProteinIdentification(List)`, `SimpleSearchEngineAlgorithm`, `ProSEAlgorithm`, `annotateSpectrum()` |
| Parameters | `Param`, and on every algorithm `getParameters()`/`setParameters()`, `getParameterValues()`/`setParameterValues({...})` |
| Misc | `getVersion()`, `getOpenMSDataPath()`, `setLogEnabled(channel, enabled)`, `getThreadingInfo()`, `getMemoryInfo()`, the file system `FS` |

`openms.d.ts` declares all of them for TypeScript, including the shapes of the plain objects
they return (e.g. `getPeaks(): { mz: Float64Array; intensity: Float32Array }`).

## Using the module

### Loading

`openms.mjs` is an ES module whose default export creates an instance of the
library. It loads `openms.wasm` from the directory it was loaded from; the file includes
OpenMS' shared data (modifications, enzymes, controlled vocabularies, ...), which the
library finds at `/share/OpenMS` of its file system. Options of
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
processed in a Web Worker, so that the page stays responsive; a worker can also read the
user's file directly, without copying it into memory first:

```js
// worker.mjs, started with new Worker("worker.mjs", { type: "module" })
import createOpenMS from "@openms/openms-wasm";
const ready = createOpenMS();                    // set onmessage before awaiting anything

onmessage = async ({ data: file }) => {          // a File from <input type="file">
  const OpenMS = await ready;
  OpenMS.FS.mkdir("/input");
  OpenMS.FS.mount(OpenMS.FS.filesystems.WORKERFS, { files: [file] }, "/input");
  const handler = new OpenMS.FileHandler();
  const exp = new OpenMS.MSExperiment();
  try {
    handler.loadExperiment(`/input/${file.name}`, exp);
    postMessage(exp.getSpectrumSummary());       // plain objects: no cleanup, cloneable
  } finally {
    handler.delete();
    exp.delete();
    OpenMS.FS.unmount("/input");
    OpenMS.FS.rmdir("/input");
  }
};
```

`examples/node/summary.mjs` shows the same with Node.js (`NODEFS`).

### Files

OpenMS reads and writes files by path. In WebAssembly, paths refer to Emscripten's
[virtual file system](https://emscripten.org/docs/api_reference/Filesystem-API.html),
`OpenMS.FS`, which keeps files in memory:

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
`storeBuffer(exp)`). Without copying the file into memory first:

- Node.js: mount a directory, `FS.mount(FS.filesystems.NODEFS, { root: "/data" }, "/data")`.
- Web Workers: mount `File` objects,
  `FS.mount(FS.filesystems.WORKERFS, { files: [file] }, "/input")`; OpenMS then reads
  them in pieces. In the threaded build, not for files that OpenMS reads on a helper thread
  while the calling thread computes (ProSE reads the spectra that way): the file system lives
  on the calling thread, so each read of the helper waits for it, and WORKERFS reads are slow.
  A ProSE search of a 78 MB mzML file on 16 threads took 29 s from WORKERFS, 2 s from a copy
  in memory (`FS.writeFile`).

Gzip- and bzip2-compressed mzML files are read as usual.

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
picked-protein FDR. `search(..., { rescore: true })` rescores the PSMs with Percolator first,
in-process as the tool's `-rescore` does (needs decoys and at least 100 PSMs; the FDR filters
then use its q-values). It fills the two lists and returns a summary: the statistics of the run
and of the fragment index, Percolator's outcome, and for open searches the modification analysis
(delta-mass histogram, mapped PTMs).

```js
const prose = new OpenMS.ProSEAlgorithm();
prose.setParameterValues({ "fragment:mass_tolerance": 0.02, "modifications:fixed": ["Carbamidomethyl (C)"] });
prose.setProgressCallback((p) => postMessage(p));   // { event, label, depth, begin, end, value }
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

`examples/prose/` is a web application around it: settings generated from ProSE's parameters,
progress, summary, score and mass-error distributions, delta masses of open searches, PSM and
protein tables, and for each PSM the sequence with its fragment coverage, the annotated spectrum
mirrored against the theoretical one and the fragment mass errors; idXML and TSV downloads. It
runs OpenMS in a Web Worker with all cores on cross-origin isolated pages. `examples/node/prose.mjs`
runs a search from the command line and reports time and memory.

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
4 GB; `low_memory` does the same for non-specific searches. A 64-bit module (Memory64) would
raise the limit (Chrome and Firefox allow up to 16 GB; Safari does not support Memory64 yet,
as far as known). Measured with a local wasm64 build, it did not pay off for ProSE: 10 to 20%
slower, 1.5 to 3 times the memory, and the semi-tryptic human search above still needed
chunking (it needs more than 16 GB without); Embind would also pass `size_t` as BigInt. More threads than about
8 to 16 hardly shortened these searches: decoy generation, parts of the index build and reading
the spectra are serial.

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

### Limits

The module is a 32-bit (wasm32) program: it can address up to 4 GB. It cannot start
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
