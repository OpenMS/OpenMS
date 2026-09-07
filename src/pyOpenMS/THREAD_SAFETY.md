# pyOpenMS Thread-Safety Contract

## Scope

pyOpenMS releases Python's Global Interpreter Lock (GIL) around selected long-running C++ calls so that other Python threads can run while the call is in progress. Releasing the GIL is an interoperability and responsiveness feature; it is **not** a general thread-safety guarantee for OpenMS objects.

This document defines the guarantee made by the pyOpenMS binding layer. General thread safety of the OpenMS core is outside its scope.

## Rules for callers

Unless an individual OpenMS API explicitly documents stronger guarantees, follow these rules:

- Do not call pyOpenMS methods concurrently on the same object unless the type or API documents concurrent read-only access. The object must be fully initialized, and no thread may write to it or receive a writable alias into it.
- Do not read or mutate an input, output, or object reachable from either concurrently with a GIL-releasing call that can modify it or expose writable aliased storage.
- Give each worker its own algorithm and data objects. Use copies or explicit synchronization when data must cross worker boundaries.
- Python callbacks routed through the streaming consumer caster (used by APIs like ImzMLFile.load with a Python consumer object) reacquire the GIL before entering Python. Not all streaming-style APIs use this caster; check the specific binding's docstring.

In particular, an Arrow table or Python buffer that aliases OpenMS storage is shared mutable state for its entire lifetime. It must not be read or written concurrently with an operation that can modify that storage.

## GIL-releasing binding inventory

The current bindings contain 62 `nb::gil_scoped_release` sites. This is an implementation inventory, not an endorsement of calling the listed operations concurrently.

| Binding source | Calls | Covered operations |
| --- | ---: | --- |
| `bindings/arrow_zerocopy.cpp` | 16 | Arrow import/export, including experimental zero-copy paths |
| `bindings/bind_experiment.cpp` | 2 | `MSExperiment.sortSpectra()` and `sortChromatograms()` |
| `bindings/bind_format.cpp` | 13 | FileHandler, imzML, and indexed mzML load/store/validation operations |
| `bindings/bind_misc.cpp` | 31 | File I/O and long-running feature-finding, deconvolution, indexing, and search operations |

**Note:** one `bind_misc.cpp` site (`TransitionListEvidenceFilter.filter`) additionally parallelizes internally via a `threads` argument; the active SWATH maps and their spectra, plus derived `candidates`/`precursor_index`, are read by multiple worker threads during the threaded scan. `transition_exp` is read before that phase to build `candidates`.

When adding a GIL-releasing binding, document any API-specific concurrency guarantee in the binding's docstring and update this inventory. Do not infer a thread-safety guarantee merely because the GIL is released.

## Avoiding OpenMP oversubscription
- Many OpenMS algorithms are internally OpenMP-parallelized, invisible from Python.
- Combining that with your own Python-level thread pool multiplies thread count (N Python workers × M OpenMP threads per call).
- Check current cap: `OpenMSBuildInfo.getOpenMPMaxNumThreads()`.
- If you need to cap oversubscription, call `OpenMSBuildInfo.setOpenMPNumThreads(n)` before entering OpenMS-parallel work. This is a task-scoped override for the current worker/context, not a single worker-global setting that automatically applies everywhere.
- After unsetting `OMP_NUM_THREADS`, each Python worker that will do OpenMS-parallel work should call `setOpenMPNumThreads(n)` itself; setting it only in the parent process does not propagate independently to workers.
- Rule of thumb: If running P Python workers concurrently, set OpenMP threads to roughly `total_cores // P` (or 1) in each worker before it starts OpenMS work. This only matters when you're doing your own Python-level threading or multiprocessing; a single-threaded script calling OpenMS functions has nothing to worry about and doesn't need to touch `setOpenMPNumThreads()`.