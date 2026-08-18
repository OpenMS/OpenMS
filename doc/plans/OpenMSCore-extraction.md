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

---

# Part 2: upstream changes that make Part 1 simpler and cleaner

Everything hard in Part 1 is hard because OpenMS keeps configuration, global state, or type
knowledge one layer below where it belongs. Each of those is independently fixable — mostly by
deleting code — and each fix is a small, separately green PR. Landing S1–S4 *before* the
extraction turns Part 1 into nearly the pure "`git mv` + build system" move the original sketch
imagined. All statements below are measured on the #9919 branch.

## S1. Make the bottom of OpenMS configuration-free

Measured: across the 18 files of §3.2, the only configure-time macros actually consumed are
`OPENMS_WINDOWSPLATFORM` (≙ `_WIN32`), `OPENMS_HAS_UNISTD_H`/`OPENMS_HAS_SYS_RESOURCE_H`/
`OPENMS_HAS_KILL` (≙ `__has_include`, C++17), `OPENMS_ASSERTIONS` (≙ `!defined(NDEBUG)` — which
is also *more* correct than today's `CMAKE_BUILD_TYPE STREQUAL Debug`, which never fires for the
Debug config of multi-config generators), and `OPENMS_NO_FLOAT_FROM_CHARS` (a compiler
feature test). `OPENMS_PRETTY_FUNCTION` is already compiler-derived. Elsewhere,
`OPENMS_BIG_ENDIAN` (4 users) is `std::endian` since C++20 and `OPENMS_HAS_PROCESS_H` has zero
users. **The core needs nothing from `configure_file` at all.**

Therefore:

* Replace the configure-time checks with compiler facts and make the platform slice of `config.h`
  a **static, checked-in header**. The §3.4 problem ("who generates `config.h` first") ceases to
  exist — there is nothing to generate.
* Replace the per-target `generate_export_header()` calls with **one small checked-in header**
  defining all library export macros (`OPENMS_DLLAPI`, `OPENMSCORE_DLLAPI`,
  `OPENSWATHALGO_DLLAPI`, `OPENMS_GUI_DLLAPI`), each keyed on the `<Target>_EXPORTS` define CMake
  already sets. One greppable file; headers work in IDEs without configuring; a new library is a
  three-line diff. This deletes §3.4 entirely, including the shared-build-include-dir plumbing
  and the generated-header special case in `install_headers()`.
* Put the rest of `config.h` on a diet separately: it currently bakes feature/dependency soup
  (`WITH_GUI`, `ENABLE_UPDATE_CHECK`, Boost/LibSVM/GLPK version numbers, `WITH_CRAWDAD`,
  COINOR/HiGHS, install paths) into a header included by every TU of every consumer — flipping
  `WITH_GUI` or the install prefix rebuilds the world. Those belong as target-level compile
  definitions on the handful of TUs that read them. Not a prerequisite for Part 1, but it is the
  same disease, and it keeps the core from owning macros about GUI and LP solvers.

## S2. Retire the push-model `GlobalExceptionHandler`

Today every `BaseException` constructor calls `GlobalExceptionHandler::getInstance().set(...)`,
pushing file/line/function/name/what into static strings that a custom `terminate()` handler
prints later. That design (from before C++11) is why `Exception.cpp` drags
`GlobalExceptionHandler` into the core — and it is independently broken: the statics are a data
race under OpenMP, the handler is only installed once the first exception has ever been
constructed, and `terminate()` prints the *last constructed* exception even when it is unrelated
to the actual termination.

The pull model deletes all of it: the terminate handler calls `std::current_exception()`,
rethrows, catches `Exception::BaseException&`, and reads the same five fields from the in-flight
object (they are all stored members). Consequences:

* the `set()`/`setMessage()` calls disappear from `Exception.cpp` and from the four outside call
  sites (`ProForma.cpp`, `TOPPBase_defs.h`, two FEATUREFINDER headers) — all mechanical deletions;
* `GlobalExceptionHandler` shrinks to a ~40-line terminate/new-handler installer that no longer
  needs to be in the core at all (`Exception.cpp` becomes dependency-free);
* the printed diagnostics get *more* accurate, and a real latent race is gone.

## S3. Localize the PID fetch in `UniqueIdGenerator`

`UniqueIdGenerator.cpp` includes `SysInfo.h` for exactly one call, `SysInfo::getProcessId()`,
used to mix the PID into the seed. A three-line `#ifdef _WIN32 GetCurrentProcessId() #else
getpid()` in the `.cpp` removes it. `SysInfo` (335 lines of psapi/mach/proc memory
introspection, 8 real users, all library-level) then stays in libOpenMS where it belongs, and
the core loses its entire platform-API surface. While in the file: the hand-rolled
pointer-plus-`omp critical` singleton can become a function-local static (thread-safe since
C++11), deleting three critical sections.

## S4. Move type knowledge out of the bottom headers — the same fix, three places

The framework's `isRealType(DataValue)` (§3.5) is one instance of a pattern; there are two more:

* **`Types.h` specializes `writtenDigits<DataValue>`** — the bottom-most header of the codebase
  forward-declares and special-cases a DATASTRUCTURES type. Measured: no library code calls
  `writtenDigits` on a `DataValue`; only `TEST_REAL_SIMILAR` does. The specialization moves into
  the same test-side support header as the `isRealType` overloads.
* **`StringUtils.h/.cpp` declare and define `toStr`/`appendToStr` for `DataValue`/`ParamValue`**
  (§3.2's wart). Root fix instead of the workaround: declare those overloads next to the types
  (`DataValue.h`/`ParamValue.h`), where every caller already has the type included. Unqualified
  callers keep working via ADL; the handful of `StringUtils::`-qualified callers get a
  deprecation shim or a one-time mechanical edit. `StringUtils` then contains *no* named OpenMS
  type, and the "3 declarations stay `OPENMS_DLLAPI` in a core header" exception in §3.2
  disappears.
* **Optionally delete the customization point altogether**: `DataValue::operator double()` and
  `ParamValue::operator double()` are implicit, so
  `std::is_floating_point_v<T> || (std::is_class_v<T> && std::is_convertible_v<T, double>)`
  classifies today's types correctly with zero registration — and covers any future library's
  variant type for free. Semantics are unchanged for current tests (a non-numeric `DataValue`
  already throws on conversion today). The support header of §3.5 then shrinks to nothing; the
  26+1 tests still add their `DataValue.h`/`ParamValue.h` includes, which they morally owed
  anyway.

Also in this family: **delete `Internal::OpenMS_locale`** (`Types.cpp`). Its initializer saves
the locale, sets "C", restores the saved locale — net effect nothing — and the global holds the
literal `"C"`. It has one reader (`OpenSwathFeatureXMLToTSV.cpp`, which can use the literal) plus
its own relic test. With it gone, `Types.cpp` is empty and is deleted; the core loses another
static-initializer global.

## S5. Finish the migration to target-based CMake

The test projects still consume `include_directories(${OpenMS_INCLUDE_DIRECTORIES})` and friends
from `CACHE INTERNAL` variables that `openms_add_library()` maintains in parallel to the targets'
real usage requirements. Every new library means more bookkeeping in that shadow system. Deleting
it (tests just `target_link_libraries(... OpenMS)`) is what makes "add a library" a non-event —
for `OpenMSCore` and for every later slice. Same direction, lower priority: the per-directory
`sources.cmake`/`sources.cmake` pairs make every file move a 4–6-list edit; one list per library
(or `CONFIGURE_DEPENDS` globs, if acceptable) would reduce §4's step 2 to the `git mv`.

## S6. Make the shipped-library list a derived fact

Part 1's §5 checklist exists because four places hardcode the library roster (wheel workflow
components, pyOpenMS copy commands, the `__init__.py` preload order, KNIME's `foreach` lists).
Two changes retire most of it: give the shipped libraries a `$ORIGIN` self-RPATH in the wheel
(the installed tree already uses `$ORIGIN/../lib`; the wheel copies bypass it — hence the
manual, order-sensitive `ctypes` preload), and drive the copy/packaging lists from the existing
`_OPENMS_EXPORT_TARGETS` list instead of naming libraries. After that, adding a library touches
the CHANGELOG and nothing else.

## Revised sequencing and what Part 1 becomes

Land as independent PRs, in any order, before the extraction: **S2, S3, S4** (each deletes more
than it adds), and **S1** at least for the platform/export headers. Then re-run Part 1, which
shrinks to:

* extraction set: **9 headers + 4 sources** (`GlobalExceptionHandler.*`, `SysInfo.*` no longer
  move; `Types.cpp` and `PrecisionWrapper.cpp` are deleted outright), every `.cpp` of which is
  dependency-free except `StringUtils.cpp` → SIMDe (`PRIVATE`, header-only);
* §3.2's declaration wart: gone (S4);
* §3.4: gone (S1);
* §3.5: a support header with two overloads — or nothing at all with the predicate variant (S4);
* §5's checklist: CHANGELOG plus whatever of S6 has not landed yet.

The original sketch claimed "8 headers, 5 sources, build-system only". After S1–S4 that claim is
essentially true — not because the plan got smaller, but because the code was changed to match
the architecture the sketch assumed.

---

# Part 3: the Catch2 answer — free the framework without creating a library

Parts 1–2 are still complex because they answer two questions at once: *free the test framework*
and *start splitting libOpenMS*. Every piece of residual machinery — the shared library, the
export macros, the `config.h` wiring, the packaging checklist, the ABI note — is a cost of the
**second** goal. Catch2 never pays any of it because it refuses that goal's premise: it depends
on nothing, so there is no library to create, version, export, or ship. This part supersedes the
sequencing above **for the framework goal**; the `OpenMSCore` extraction remains the right first
slice *of the split*, on the split's own timeline and justification.

## 3.1 The load-bearing trick, named

Catch2's cleanliness is not "it avoids dependencies by being careful". It is one structural rule:

> **The framework's *compiled* code touches only `std` types. Everything type-specific happens at
> the macro-expansion site, inside the user's TU, where the user's headers are visible.**

Knowledge flows *into* the framework by ADL/`operator<<`/trait specialization/registration —
never by the framework naming a type. The sketch's objection ("we can't be std-only —
`TEST_EQUAL` and `TEST_REAL_SIMILAR` legitimately format OpenMS types") conflated *tests
formatting OpenMS types* (true, and it happens at the expansion site) with *the framework linking
OpenMS formatting* (not needed). Measured proof that OpenMS is already halfway there:
`testEqual`'s failure output streams `expression_1`/`expression_2` with bare `operator<<` — the
expansion-site mechanism, working today.

## 3.2 The five remaining couplings, mapped to Catch2's mechanisms

| coupling (measured) | Catch2 mechanism | OpenMS change |
|---|---|---|
| `StringUtils::toStr` in `testEqual`'s `std::string == numeric` quirk (reproduces old `String(114) == "114"`) | `StringMaker` — the framework owns its stringification | framework-local numeric formatter; **the one place needing care**: for floating point it must format like `toStr` or string-vs-double comparisons shift — the 720-test suite is the oracle (~20 lines) |
| `ListUtilsIO.h` include (so `TEST_EQUAL` can print vectors) | users provide `operator<<`; framework prints "{?}" otherwise | framework-local `operator<<`-if-streamable printer, or tests include `ListUtilsIO.h` themselves — either way the include leaves `ClassTest.h` |
| `BaseException` caught/printed by `ClassTest.cpp`; `ConversionError`/`IndexOverflow` ctors referenced via inline `StringUtils.h` code | catch `std::exception` + registered translators | `BaseException` already **is** a `std::runtime_error`; the framework catches `std::exception` only, and name/file/line context comes from a **registered exception translator** (changing `what()` itself is off the table — 29 tests compare it verbatim in `TEST_EXCEPTION_WITH_MESSAGE`; see the implementation spec). `TEST_EXCEPTION(Exception::Precondition, …)` is untouched: its typed `catch` expands in the test TU. The ctor references disappear with the `StringUtils.h` include |
| `StringUtilsHelper::extractDouble` + `trim` in `FuzzyStringComparator`/`ClassTest.cpp` | the framework owns its own parsing | framework-local copies (49 + ~6 lines, incl. the `from_chars`-for-float fallback). Parsing arbitrary file text is semantically independent of the library by nature; drift risk is confined here |
| `UniqueIdGenerator::setSeed` at `START_TEST` | event listeners — registration by the side that knows | a support TU in the openms class-test project whose static initializer calls `setSeed(2453440375)` (one `target_sources` line; runs even earlier than today's `START_TEST` call). OpenSwathAlgo tests simply don't compile it. This re-adds what #9919's hook-removal took out — but as *test-project code*, not framework API |
| `isRealType`/`writtenDigits<DataValue>` | type traits with sensible defaults | the S4 predicate (`is_floating_point` ∨ class-convertible-to-`double`) — no registration at all |

## 3.3 What this yields, and what it costs

`OpenMSTestFramework` becomes a **std-only static library**: zero OpenMS includes, zero link
dependencies, roughly 200–250 owned utility lines. `src/testframework/CMakeLists.txt` loses the
entire `$<TARGET_PROPERTY:OpenMS,…>` propagation block *and* needs no replacement; the
OpenSwathAlgo tests drop `"$<LINK_LIBRARY:needed,OpenMS>"`; the acid test holds. **No new shipped
artifact, no `OPENMSCORE_DLLAPI`, no config wiring, no packaging checklist, no ABI note.**

Honest costs:

* ~250 lines duplicated in spirit (`trim`, number extraction, numeric formatting). The failure
  mode is loud (test failures), not silent, and confined to the fuzzy comparator's parsing and
  the string-vs-numeric `TEST_EQUAL` quirk.
* `what()` gaining context changes the text of uncaught-exception reports tool-wide — an
  improvement, but a visible one.
* The framework still cannot validate XML or seed anything by itself — by design; #9919 already
  moved those to the right side of the line.

## 3.4 Revised plan of record

1. **PR: framework goes std-only** — local utils + formatter, `std::exception` catch + registered
   exception translators (a `what()`-completion PR was considered and dropped: 29 tests compare
   `what()` verbatim in `TEST_EXCEPTION_WITH_MESSAGE`), includes dropped, CMake decoupling,
   seeding support TU, OpenSwathAlgo link line cleaned. Acid test.
2. **PR: the S4 predicate** with the 26+1 test includes.

None of these touches packaging, installers, or wheels. `OpenMSCore` (Parts 1–2) proceeds — or
doesn't — purely as the first slice of the libOpenMS split, no longer coupled to the framework.

*Why not simply adopt Catch2 itself?* 720 test files speak `TEST_*`, and the fuzzy file
comparison plus CTest wiring live here. Wrapping the `TEST_*` macros over Catch2 is feasible as a
separate, larger migration; this part gets Catch2's *architecture* — which is what makes the
framework split-proof — without touching a single existing test.

**Implementer handoff:** the file-by-file specification for this part, with copy sources,
verification steps, and acceptance criteria, is
[`testframework-std-only-implementation.md`](testframework-std-only-implementation.md) — based on
`develop`, so it does not require #9919 to merge first (#9919's branch serves as a donor for
three verified hunks). Details discovered while writing it revise §3.4 above: the `isRealType`
predicate must be part of the main PR (the header cannot go std-only while it still includes
`DataValue.h`); the `what()`-completion PR is dropped — 29 test files compare `what()` verbatim
in `TEST_EXCEPTION_WITH_MESSAGE`, so exception naming is injected via a registered translator
(the exact Catch2 `CATCH_TRANSLATE_EXCEPTION` design) from the same support TU that seeds the ID
generator; and sequencing std-only *before* extraction makes #9919's `needed` link feature,
`LINK_LIBRARY_OVERRIDE`, CMake 3.24 bump, and OpenMS→framework usage-requirement propagation all
unnecessary.
