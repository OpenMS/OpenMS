# Plan: making `OpenMSTestFramework` independent of libOpenMS

Follow-up to [#9919](https://github.com/OpenMS/OpenMS/pull/9919), which extracted the class-test
framework into the static library `OpenMSTestFramework` but left it using a handful of libOpenMS'
low-level utilities. This document turns the sketch from the PR discussion into an implementation
plan, verified against the code on the #9919 branch. Everything quantified below was measured
there, not estimated.

**Goal (acid test):** `Datastructures_test` — an OpenSwathAlgo class test — links and passes with
no `OpenMS` on its link line, and
`src/tests/class_tests/openswathalgo/CMakeLists.txt` drops its explicit
`"$<LINK_LIBRARY:needed,OpenMS>"`.

**Mechanism:** extract the utilities the framework uses into a new shared library `OpenMSCore` —
the first real slice of splitting libOpenMS — and replace the framework's only knowledge of
higher-layer types (`DataValue`/`ParamValue` in `TEST_REAL_SIMILAR`) with a customization point on
the test side. The framework then depends on `OpenMSCore` only, and so does every future split-off
library, which is what makes the framework "work with split libraries".

---

## 1. What the framework actually needs (measured)

Compiling `ClassTest.cpp` and `FuzzyStringComparator.cpp` standalone (only the framework and
libOpenMS include trees plus stand-ins for the two generated headers) succeeds, and `nm -uC` on the
two objects gives the complete link-level surface the future core must export for the framework:

| area | undefined symbols |
|---|---|
| `Exception` | `BaseException` dtor, `getFile/getFunction/getLine/getName`, **`typeinfo for BaseException`**, `ConversionError` ctor, `IndexOverflow` ctor |
| `StringUtils` | `appendToStr(int, std::string&)`, `StringUtilsHelper::extractDouble` |
| `UniqueIdGenerator` | `setSeed(UInt64)` |

(Everything else from `StringUtils.h`/`ListUtils.h`/`PrecisionWrapper.h`/`Types.h` the framework
uses is header-inline.) The `typeinfo` entry matters: `TEST_EXCEPTION` catches
`Exception::BaseException` thrown *inside libOpenMS* from a handler *in the test executable*, so
the exception's RTTI must be a single shared definition. That constrains the design (§3.1).

## 2. Corrections to the sketch in the #9919 discussion

The PR-comment sketch measured the *header* closure (8 headers, 5 sources, ~1,150 lines, "no
source edits, build-system only"). The **implementation** closure is larger, and the extraction is
not source-edit-free:

| file the sketch listed | what its `.cpp` additionally drags in |
|---|---|
| `Exception.cpp` | `CONCEPT/GlobalExceptionHandler.h` (every exception ctor calls `GlobalExceptionHandler::getInstance().set…()`) |
| `UniqueIdGenerator.cpp` | `SYSTEM/SysInfo.h` (`SysInfo::getProcessId()` is mixed into the seed) |
| `StringUtils.cpp` | `DATASTRUCTURES/DataValue.h`, `DATASTRUCTURES/ParamValue.h` (three `toStr`/`appendToStr` overload definitions), `SYSTEM/SIMDe.h` (SIMD `trim` implementation — a vendored third-party header) |
| `FuzzyStringComparator.cpp` | `SYSTEM/PathUtils.h` (header-only, already known) |

Consequences, all handled below: `GlobalExceptionHandler` and `SysInfo` join the core;
the three `DataValue`/`ParamValue` overload *definitions* move up into libOpenMS (§3.2); the core
compiles against SIMDe (`PRIVATE`, header-only — the *framework* still needs no third-party
headers, that claim was about the framework and remains true).

Confirmed as stated in the sketch:

* The trait cost is exactly **26 test files** (23 relying on `ClassTest.h` for `DataValue.h`,
  3 for `ParamValue.h`), plus 1 in `openms_gui`.
* Nothing in the closure needs Boost or Qt (`config.h` defines `OPENMS_PRETTY_FUNCTION` with
  compiler builtins; no core `.cpp` uses `Macros.h`).
* No seeding hook and no exception translator are needed: `UniqueIdGenerator` and `Exception` are
  both in the core, so the framework keeps calling/catching them directly.

## 3. Design decisions

### 3.1 `OpenMSCore` is a **shared** library

A static core linked into both libOpenMS and (via the framework) each test executable would give a
test process two copies of every core global. Four of those are load-bearing:

1. **`UniqueIdGenerator::instance_`** — `START_TEST` seeds it; libOpenMS code draws unique IDs
   from it. With two copies, the framework seeds the executable's instance while the library uses
   its own unseeded one → IDs written into test output files become nondeterministic → fuzzy
   reference mismatches. On Windows there is no symbol interposition, so this happens always; on
   Linux, `openms_hide_static_archive_symbols()` from #9919 (`--exclude-libs,ALL`) would hide the
   executable's copy from the dynamic table and break the interposition that would otherwise paper
   over it.
2. **`GlobalExceptionHandler`'s singleton** — same duplication.
3. **`Internal::OpenMS_locale` in `Types.cpp`** — a static initializer that sets up the "C" locale
   handling (pyOpenMS depends on its exact behavior); it would run once per module.
4. **`Exception::BaseException` vtable/typeinfo** — duplicated RTTI breaks cross-library `catch`
   under `-fvisibility=hidden` on macOS (typeinfo compared by pointer).

This is the Thrift double-destruction lesson from #9919 generalized, so: shared library, one copy
per process. libOpenMS links it `PUBLIC` (its own headers expose `Exception`, `StringUtils`, …).
The cost is a new runtime artifact, which is a packaging checklist (§5), not an architecture
problem — `openms_add_library()` already automates install, export, and Windows DLL copying, and
`OpenSwathAlgo` is the shipping precedent.

### 3.2 Scope: the minimal closure; `DataValue`/`ParamValue` stay out

Files that move (`git mv`, history follows), 11 headers + 7 sources, ~4,400 lines:

| layer | headers | sources |
|---|---|---|
| `CONCEPT` | `Types.h`, `Exception.h`, `GlobalExceptionHandler.h`, `PrecisionWrapper.h`, `UniqueIdGenerator.h` | `Types.cpp`, `Exception.cpp`, `GlobalExceptionHandler.cpp`, `UniqueIdGenerator.cpp` |
| `DATASTRUCTURES` | `StringUtils.h`, `ListUtils.h`, `ListUtilsIO.h`, `TypeAliases.h` | `StringUtils.cpp`, `ListUtils.cpp` |
| `SYSTEM` | `SysInfo.h`, `PathUtils.h` | `SysInfo.cpp` |

(`PrecisionWrapper.cpp` is an empty namespace — delete it instead of moving it.) All install
paths stay `<prefix>/include/OpenMS/...`, so **no consumer include changes anywhere**; the
`install_headers()` machinery already merges per-library include trees.

`DataValue`/`ParamValue` measure as core-eligible themselves (their only extra dependency is the
std-only `HashUtils.h`; +~2,800 lines), but they stay out of this slice: it keeps the first slice
small, and it preserves the point of the exercise — the *framework* must not know concrete library
value types, so that `TEST_REAL_SIMILAR` works for any future library via the customization point
rather than by enumeration. They can join the core in a later slice without rework (see below).

The one genuine source consequence: `StringUtils.cpp` currently defines
`appendToStr(const DataValue&, bool, std::string&)`, `toStr(const DataValue&, bool)` and
`toStr(const ParamValue&, bool)` (~45 lines). These move to `DataValue.cpp` / `ParamValue.cpp` in
libOpenMS. Their declarations (and the inline `StringConversions` forwarders) stay in
`StringUtils.h` against the existing forward declarations, but keep the `OPENMS_DLLAPI` macro —
the only OpenMS-DLL symbols declared in a core header, with a comment saying so. This is
deliberately self-healing: when `DataValue`/`ParamValue` later join the core, the definitions ride
along in their `.cpp`s and only the macro flips.

### 3.3 Export macro: `OPENMSCORE_DLLAPI` via the existing machinery

`openms_add_library(TARGET_NAME OpenMSCore … DLL_EXPORT_PATH "OpenMS/")` already generates
`OpenMS/OpenMSCoreConfig.h` defining `OPENMSCORE_DLLAPI` (`generate_export_header`, same as
`OPENSWATHALGO_DLLAPI`). Keeping the name `OPENMS_DLLAPI` inside the moved headers is not an
option: a symbol's macro must expand to dllexport exactly when its *own* library is being built,
and libOpenMS keeps `OPENMS_DLLAPI` for its remaining 700+ headers.

Mechanical edit in the 11 moved headers: 77 `OPENMS_DLLAPI` occurrences → 74 become
`OPENMSCORE_DLLAPI` (3 stay, §3.2), and their `#include <OpenMS/OpenMSConfig.h>` lines become
`<OpenMS/OpenMSCoreConfig.h>`.

### 3.4 `config.h` ownership and the `OpenMSConfig.h` transitive include

Constraint, measured: 730 libOpenMS headers use `OPENMS_DLLAPI` but only 54 include
`OpenMSConfig.h` directly — the rest get it via `config.h` line 25. So `config.h` must keep
providing `OPENMS_DLLAPI`, and `OpenMSConfig.h` must keep existing.

Recommended wiring:

* `configh.cmake` + `config.h.in` move to the core subproject (`config.h` is included by the
  core's own `Types.h`; the core owns platform detection). `config.h` keeps its
  `#include <OpenMS/OpenMSConfig.h>`.
* The generated per-target export headers move to one shared build include directory
  (`generate_export_header(... EXPORT_FILE_NAME ${OPENMS_HOST_BINARY_DIRECTORY}/include/...)`) on
  every library's `PUBLIC` `BUILD_INTERFACE`, so the core (and the framework) see
  `OpenMSConfig.h` without referencing the openms subproject's binary dir. Smallest-diff
  alternative: the core adds `src/openms`'s binary include dir to its build interface — functional
  because all generation happens at configure time, before anything compiles (the same property
  #9919's `$<TARGET_PROPERTY:...>` propagation relies on today), just less tidy.
* Installed trees are unaffected either way: both generated headers install into
  `include/OpenMS/` as before.

This is the fiddliest part of step 1; it is pure build-system plumbing.

### 3.5 The customization point for `TEST_REAL_SIMILAR`

`ClassTest.h` currently includes `DataValue.h` solely for two overloads
(`isRealType(const DataValue&)` / `isRealType(const ParamValue&)`). Because the macros call
`TEST::isRealType(...)` *qualified*, the extension mechanism must live in
`OpenMS::Internal::ClassTest` — either extra overloads or a trait the macro switches to; extra
overloads are the smaller diff. Following the `TestFileValidation.h` precedent from #9919, a new
test-side header in the openms class-test project (e.g.
`src/tests/class_tests/openms/include/OpenMS/ClassTestRealTypes.h`) includes
`DataValue.h`/`ParamValue.h` and provides the two overloads verbatim; `ClassTest.h` drops the
include and the overloads (and its now-vestigial `<OpenMS/OpenMSConfig.h>` include — the framework
exports nothing).

* Measured cost: the 26 openms tests (+1 GUI test) that today get `DataValue.h`/`ParamValue.h`
  transitively add one include.
* Fallout containment: dropping `DataValue.h` from `ClassTest.h` would also silently drop
  `TypeAliases.h` (108 openms tests mention `StringList`/`IntList`/`DoubleList`, mostly satisfied
  via their class-under-test headers — but not provably all). Since `TypeAliases.h` is a core
  header, `ClassTest.h` simply keeps including it explicitly (it is also the honest statement:
  `TEST_EQUAL` prints those lists via `ListUtilsIO.h`). That pins the expected fallout at the
  measured 27 files; anything beyond that surfaces as a trivial missing-include compile error.

## 4. Steps

Two PRs; each independently green. (A third, optional, later.)

### PR 1 — extract `OpenMSCore` (no consumer-visible change)

1. New subproject `src/core/` (added before `openswathalgo` in `src/CMakeLists.txt`), built with
   `openms_add_library(TARGET_NAME OpenMSCore … DLL_EXPORT_PATH "OpenMS/")`;
   `PRIVATE` deps: `SIMDe`, `OpenMP::OpenMP_CXX` if enabled (for `UniqueIdGenerator`'s
   `omp critical`); no public third-party deps. Call `openms_doc_path()` like the other libs.
   (Target names in `target_link_libraries` resolve at generate time, so the `SIMDe` target coming
   from `src/openms/extern` still works; hoist that `extern` handling only if it bothers anyone.)
2. `git mv` the files of §3.2; delete their entries from
   `src/openms/source/{CONCEPT,DATASTRUCTURES,SYSTEM}/sources.cmake` and
   `src/openms/include/OpenMS/{CONCEPT,DATASTRUCTURES,SYSTEM}/sources.cmake`.
3. Macro rename + export-header include swap in the moved headers (§3.3), `config.h` wiring
   (§3.4), move the three `DataValue`/`ParamValue` overload definitions (§3.2).
4. `target_link_libraries(OpenMS PUBLIC OpenMSCore)`. `OpenSwathAlgo` stays independent — it must
   not grow a core dependency.
5. Packaging checklist (§5).

Verify: full CI matrix; `nm -D libOpenMS.so` no longer defines `Exception`/`StringUtils`/
`UniqueIdGenerator` symbols and `readelf -d` shows `NEEDED libOpenMSCore.so`; wheels build and the
20 wheel-install jobs pass; `UniqueIdGenerator_test` and one ID-writing format test confirm
seeding still reaches the library (the §3.1 hazard, ruled out by the shared design but worth the
explicit check on Windows CI).

### PR 2 — cut the framework's cord

1. `ClassTestRealTypes.h` + the 27 includes; `ClassTest.h` drops `DataValue.h` and the two
   overloads, keeps `TypeAliases.h` (§3.5).
2. `src/testframework/CMakeLists.txt`: replace the `$<TARGET_PROPERTY:OpenMS,...>` include/defs
   propagation and `add_dependencies(OpenMSTestFramework OpenMS)` with
   `target_link_libraries(OpenMSTestFramework PUBLIC OpenMSCore)`. The #9919 Thrift/Arrow hazard
   does not reappear: `libOpenMSCore.so` contains no Arrow symbols, and the openms tests' explicit
   ordering of `OpenMS` before the static Arrow archives is untouched. No `needed` feature for the
   core either: `ClassTest.o` is pulled into every test by `START_TEST` and references core
   symbols, so `--as-needed` keeps the core by construction.
3. `src/tests/class_tests/openswathalgo/CMakeLists.txt`:
   `target_link_libraries(${i} OpenSwathAlgo OpenMSTestFramework)` — the explicit
   `"$<LINK_LIBRARY:needed,OpenMS>"` and its comment go away. The openms/openms_gui test link
   lines (and `LINK_LIBRARY_OVERRIDE_OpenMS`) stay as they are: those tests genuinely test
   libOpenMS.
4. Acid test: the OpenSwathAlgo tests link with no `OpenMS` anywhere in their link closure and
   pass. Optionally assert it permanently with a one-line CTest
   (`readelf -d Datastructures_test | grep -c libOpenMS.so` = 0 on ELF).

### Later slices (explicitly out of scope now)

* `DataValue`/`ParamValue`/`HashUtils` into the core (the three overload definitions come home;
  `OPENMS_DLLAPI` → `OPENMSCORE_DLLAPI` on them). Do this when the next split slice wants it, not
  before.
* Further libOpenMS splitting; the core's dependency list is its seed specification.
* A fully std-only framework (Catch2 purism) stays rejected: it would mean duplicating
  `toStr`/`extractDouble`/trim behavior the reference-file comparisons depend on. The inversion
  goal is met at the knowledge level: bottom-layer utilities plus one test-side customization
  point.

## 5. Ripple checklist for the new shared library

Automatic via `openms_add_library()` / `install_library()` — no edits: install +
`OpenMSTargets.cmake` export (build tree and installed), `RUNTIME_DEPENDENCY_SET`, Windows DLL
copying to test/doc bins, `$<TARGET_RUNTIME_DLLS>` copying, doxygen path. CMake consumers
(`find_package(OpenMS)`, TOPP tools, GUI, pyOpenMS link lines, `FuzzyDiff`) get the core
transitively through `OpenMS`'s `PUBLIC` link.

Needs explicit edits:

| file | edit |
|---|---|
| `.github/workflows/pyopenms-wheels-cibuildwheel.yml` | add `cmake --install . --component OpenMSCore_headers` next to the existing `OpenMS_headers`/`OpenSwathAlgo_headers` lines (4 platforms; the `ninja OpenMS OpenSwathAlgo` lines are fine — the core builds as a dependency, and `--component library` installs it) |
| `src/pyOpenMS/CMakeLists.txt` | third target variable next to `OPENMS_TARGET`/`OPENSWATHALGO_TARGET` (~l. 74); add `$<TARGET_FILE:OpenMSCore>` to the `pyopenms_copy_deps` commands (~l. 412) and to `_pyopenms_dll_dirs` (~l. 539) |
| `src/pyOpenMS/pyopenms/__init__.py` | Linux wheel preload: `libOpenMSCore.so` **before** `libOpenSwathAlgo.so`/`libOpenMS.so` (~l. 170) |
| `cmake/knime_package_support.cmake` | add `OpenMSCore` to the two `foreach(... OpenMS OpenSwathAlgo)` lists (l. 326, 368) |
| `CHANGELOG` | note the new library and the link-level break |
| `AGENTS.md` / docs | mention the library layout where the existing libs are listed |

macOS `mac_fix_dependencies.rb` and delocate, and Linux auditwheel, operate on discovered dylibs
generically — verify in the wheel jobs, expect no edits.

## 6. Cost and compatibility notes

* **Compile time** (asked in the PR thread): net TU count −1 (`PrecisionWrapper.cpp` deleted);
  everything else is the same code compiled once, in a different target. Rebuild granularity
  improves at the link level (touching core internals relinks core, not all of libOpenMS —
  header edits still fan out exactly as before). Test binaries are unchanged; the new
  `libOpenMSCore` is a few hundred KB.
* **Runtime**: none. Same code, one call's worth of PLT indirection on cross-library calls.
* **ABI / downstream**: link-level breaking — the moved symbols leave `libOpenMS.so`. CMake
  consumers are unaffected (transitive `PUBLIC` link); anyone linking `-lOpenMS` by hand adds
  `-lOpenMSCore`; distro/conda packaging picks up one more library. Between minor releases, and
  in a cycle that already removed `OpenMS::String` from the ABI, this is acceptable — CHANGELOG
  entry regardless.
* **Behavior**: none intended anywhere — no formatting, comparison, validation, or seeding
  semantics change. The trait/overload move is observable only as a compile error for a test that
  uses `TEST_REAL_SIMILAR` on a new variant-like type without providing the overload, which is the
  designed extension point.
