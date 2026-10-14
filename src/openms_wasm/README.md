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
| Identification | `PeptideHit`, `PeptideIdentification(List)`, `ProteinHit`, `ProteinIdentification(List)`, `SimpleSearchEngineAlgorithm` |
| Parameters | `Param`, and on every algorithm `getParameters()`/`setParameters()`, `getParameterValues()`/`setParameterValues({...})` |
| Misc | `getVersion()`, `getOpenMSDataPath()`, `setLogEnabled(channel, enabled)`, the file system `FS` |

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
  them in pieces.

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

### Errors

An OpenMS exception becomes a JavaScript `Error` whose message names it, e.g.
`ParseError: ... while parsing 'PEPT[X'`.

### Threads

`@openms/openms-wasm` runs on a single thread: the parts of OpenMS parallelized with OpenMP
run sequentially. `@openms/openms-wasm-mt` is the same module (same files, same API) built
with threads: OpenMP regions (e.g. mzML decoding, feature finding, database search) and the
helper threads of OpenMS run on WebAssembly threads, Web Workers (Node.js: worker threads)
that share the module's memory.

```js
import createOpenMS from "@openms/openms-wasm-mt";
const OpenMS = await createOpenMS({ numThreads: 4 });   // optional
OpenMS.getThreadingInfo();   // { threaded: true, available: true, numThreads: 4 }
```

- **Node.js** needs nothing else. Without `numThreads`, OpenMS uses `OMP_NUM_THREADS` of the
  process, else all logical cores.
- **Browsers** allow shared memory only on cross-origin isolated pages: serve the page with
  `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`
  (`credentialless` also works in Chrome and Firefox, not in Safari). Elsewhere the module
  runs single-threaded (`available: false`, `numThreads: 1`). The default there is half
  the logical cores, at most 4; the Web Workers for these threads are started with the
  module. Load the module in a Web Worker: on the page's main thread, waiting threads
  busy-wait and block the page. Browser use of the threaded module is not validated yet.

`getThreadingInfo()` exists in both packages (`threaded: false` in `@openms/openms-wasm`).

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

The multi-threaded package (`@openms/openms-wasm-mt`) comes from the preset
`wasm32-emscripten-mt-release` (build directory `build/wasm32-emscripten-mt-release`): it sets
`OPENMS_WASM_THREADS=ON` and builds the dependencies a second time, for the
`wasm32-emscripten-mt` triplet (with `-pthread`), as a module with threads can only be linked
from objects compiled for it.

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
them pass except those that exercise what the module does not have: threads
(`GlobalExceptionHandler`, `LogConfigHandler`, `LogStream`, `UniqueIdGenerator`, `File`,
`TOPPBase`: they start `std::thread`s), other programs and the network (`ExternalProcess`,
`TOPPExternalToolBase`, `Network`; `File` also looks for programs next to its own), and the
TOPP tool registry (`ToolHandler`: it lists the tools of the build, and the WebAssembly
build has none).

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
