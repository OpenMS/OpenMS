# pyOpenMS Thread-Safety Contract

## Scope

pyOpenMS releases Python's Global Interpreter Lock (GIL) around selected long-running C++ calls so that other Python threads can run while the call is in progress. Releasing the GIL is an interoperability and responsiveness feature; it is **not** a general thread-safety guarantee for OpenMS objects.

This document defines the guarantee made by the pyOpenMS binding layer. General thread safety of the OpenMS core is outside its scope.

## Rules for callers

Unless an individual OpenMS API explicitly documents stronger guarantees, follow these rules:

- Do not call pyOpenMS methods concurrently on the same object.
- Do not read or mutate an input, output, or object reachable from either while a GIL-releasing call is using it. This includes objects obtained through reference-returning accessors and writable NumPy, memoryview, or Arrow views.
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
