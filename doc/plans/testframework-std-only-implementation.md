# Implementation spec: std-only `OpenMSTestFramework`, based on `develop`

Handoff document for an implementer (human or agent). **Base: `develop`** (verified at
`c77f588`), where the class-test framework still lives inside libOpenMS
(`src/openms/{include/OpenMS,source}/CONCEPT/`: `ClassTest.h/.cpp`,
`FuzzyStringComparator.h/.cpp`, `MacrosTest.h`). This program does **not** require PR
[#9919](https://github.com/OpenMS/OpenMS/pull/9919) to merge; it reaches the same end state by a
cleaner route (see "Relationship to #9919" at the end). The #9919 branch
(`origin/claude/openms-test-framework-refactor-ssutqy`) is used as a **donor** for three verified
hunks — fetch it, but do not base on it. All file:line references below are to `develop`.
Background and rationale: Part 3 of [`OpenMSCore-extraction.md`](OpenMSCore-extraction.md).

**Objective.** The class-test framework becomes its own static library `OpenMSTestFramework`
(`src/testframework/`) that compiles and links against the C++ standard library only. OpenMS
behavior (exception naming, unique-ID seeding) is injected from the test projects via
registration.

**Acceptance (in order):**
1. After PR 1: `src/openms/source/CONCEPT/ClassTest.cpp` and `FuzzyStringComparator.cpp` have no
   `#include <OpenMS/...>` other than the framework's own five headers, and compile unchanged
   inside libOpenMS.
2. After PR 2: `nm -uC libOpenMSTestFramework.a` shows **no** undefined `OpenMS::` symbols.
3. After PR 2: the five OpenSwathAlgo tests link as `target_link_libraries(${i} OpenSwathAlgo
   OpenMSTestFramework)` — no `OpenMS` — and pass (`readelf -d Datastructures_test` shows no
   `libOpenMS.so`).
4. Full `ctest` green after every PR, including the ≈1600 `TOPP_FuzzyDiff`-based TOPP tests (the
   oracle for all copied parsing/formatting code) — with **no** change to any reference file.
5. Failure output stays human-equivalent: values print; unexpected OpenMS exceptions in openms
   tests still report their name/file/line (translator path).

**Ground rules — do NOT:**
* change `BaseException::what()`/`getMessage()` (29 test files compare `what()` verbatim in
  `TEST_EXCEPTION_WITH_MESSAGE`);
* change library formatting/parsing (`StringUtils`, `Types.h` — except the PR 3 deletion);
* bump `cmake_minimum_required` (stays 3.21, `CMakeLists.txt:10`) or introduce link features
  (`$<LINK_LIBRARY:...>`) — this plan is designed so they are unnecessary;
* rename installed header paths (`<OpenMS/CONCEPT/ClassTest.h>` must keep working for ~720 tests);
* touch reference files or TOPP tools (other than `FuzzyDiff`'s link line in PR 2).

Three PRs, each independently green.

---

## PR 0 — make tmp-file validation explicit (donor: #9919)

On `develop`, `END_TEST` automatically schema-validates `NEW_TMP_FILE` outputs: `ClassTest.cpp`
includes ten `FORMAT/*` headers (`:16–:25`), implements `validate()` (`:243`, using
`FileHandler::getType` and per-format `isValid`), and calls it from `endTestPostProcess`
(`:730`). A std-only framework cannot contain this, and #9919 already solved it by measurement —
take its hunks verbatim:

1. `src/tests/class_tests/openms/include/OpenMS/TestFileValidation.h` (header-only helper,
   `VALIDATE_TMP_FILES` / `VALIDATE_FILE`) — copy from the donor branch.
2. The `VALIDATE_TMP_FILES` call added to the **26 tests that actually validate files** — the
   donor's diff *is* the measured list (all 57 `NEW_TMP_FILE`+XML tests were run against the
   automatic-validation build; 26 validated something, each now validates exactly the same
   files). Copy those 26 test-file hunks.
3. Delete from `ClassTest.cpp`: the ten `FORMAT/*` includes, `validate()` and its call, the
   `tmp_file_list`-validation plumbing the donor deletes.

Behavioral note for the PR text (same as #9919's): validation becomes opt-in for *new* tests; the
`VALIDATE_TMP_FILES` docs say when to add it. Verify: full ctest; the 26 tests still fail if fed
a deliberately broken XML (spot-check one, e.g. a featureXML test, by corrupting its tmp output
locally).

## PR 1 — the framework goes std-only, in place

Framework files stay in `src/openms/` and stay compiled into libOpenMS (only their own two
`.cpp`s include the framework headers — verified). Tests keep linking `${OpenMS_LIBRARIES}`, so
everything stays green while the sources lose every OpenMS dependency. This is the substantive
PR; land as staged commits in this order.

### 1.1 New header `src/openms/include/OpenMS/CONCEPT/ClassTestUtils.h`

Framework-owned utilities, namespace `OpenMS::Internal::ClassTest` (+ nested `detail`). Add to
`include/OpenMS/CONCEPT/sources.cmake`. Contents:

1. **`writtenDigits<T>()`** — framework's own: floating `T` → `std::numeric_limits<T>::digits10`;
   integral → `digits10`; class convertible to `double` → `digits10` of `double`; default 6
   (matches `Types.h`'s primary template). Required because the OpenSwathAlgo tests use
   `TEST_REAL_SIMILAR` (3 files) and must not need `Types.h`.
2. **`isRealType`** — replace the five overloads (`ClassTest.h:61–:89`) with one template:
   `std::is_floating_point_v<C> || (std::is_class_v<C> && std::is_convertible_v<const C&, double>)`,
   `C = std::remove_cvref_t<T>`. `DataValue::operator double()` and `ParamValue::operator
   double()` are implicit (verified), so both classify `true` with zero registration;
   `std::string`/containers/pointers stay `false`. `TEST_REAL_SIMILAR` (`:670`) keeps calling
   `TEST::isRealType(a)` unchanged.
3. **`detail::toString(T)`** for arithmetic `T` — **copy verbatim** `appendNumeric` from
   `StringUtils.cpp:71` (NaN→`"NaN"`, ±inf→`"inf"`/`"-inf"`, `std::to_chars`
   shortest-round-trip, the `1e4`/`1e-2` scientific thresholds, keep-one-fractional-digit rule)
   plus the integer `to_chars` path. Must byte-match `StringUtils::toStr` for numeric arguments —
   the `TEST_EQUAL(std::string, numeric)` quirk depends on it. Do not improve it; comment the
   source of truth at both sites.
4. **`detail::printValue(std::ostream&, const T&)`** — if streamable → `os << v`; else if range
   of streamables → `[e1, e2, ...]` (separator `", "`, matching `ListUtilsIO.h:29`); else
   `"<unprintable>"`. A named function, **not** an `operator<<` — a framework `operator<<` for
   `std::vector<T>` would be ambiguous with ADL-found `OpenMS::operator<<` when tests include
   `ListUtilsIO.h`.
5. **`detail::trim/split/join`** — plain implementations; `split` reproduces
   `ListUtils::create<std::string>` semantics (13-line source in `ListUtils.h`; FuzzyDiff
   whitelist TOPP tests are the oracle).

### 1.2 `ClassTest.h`

* Includes (`:15–:28`): delete `PrecisionWrapper.h`, `Types.h`, `StringUtils.h`, `DataValue.h`,
  `OpenMSConfig.h`, `config.h`, `ListUtilsIO.h`; add `ClassTestUtils.h`. (`MacrosTest.h` stays.)
* Delete the five `isRealType` overloads (`:61–:89`).
* `testEqual`: `StringUtils::toStr(expression_2)` (`:289`) → `detail::toString(...)`; the bare
  `stdcout << expression_N` failure prints (in `testEqual`/`testNotEqual`) →
  `detail::printValue(stdcout, ...)`. Keep the `XMLCh*` and enum `if constexpr` branches as-is.
* `TEST_REAL_SIMILAR` (`:670`): `writtenDigits(a)` → `TEST::writtenDigits(a)` (twice).
* Exception macros: in `TEST_EXCEPTION` (`:806`) delete the
  `catch (::OpenMS::Exception::BaseException& e)` clause (`:819`) and collapse the fallbacks to
  `catch (...) { TEST::exception = 2; TEST::exception_name = TEST::describeCaughtException(); }`;
  same surgery in `TEST_EXCEPTION_WITH_MESSAGE` (`:928`, clause `:946`) — its typed
  `et.what() == message` clause is untouched. Adjust the outcome-printing switches accordingly.
* Precondition macros: `#ifdef OPENMS_ASSERTIONS` (`:888`, `:905`) → `#ifndef NDEBUG` (14 test
  files use `TEST_PRECONDITION_VIOLATED`; `!NDEBUG` tracks the test TU's real compile mode,
  which is what must agree with the library's `OPENMS_PRECONDITION`).
* `OPENMS_WINDOWSPLATFORM` → `_WIN32` wherever it appears in framework files.
* New declarations (see 1.3): `registerExceptionTranslator`, `describeCaughtException`. While the
  framework is still inside libOpenMS, decorate the new `.cpp`-backed functions `OPENMS_DLLAPI`
  like their neighbors (PR 2 strips all of these).

### 1.3 `ClassTest.cpp`

* Replace: `File::remove` (`:192`) / `File::exists` (`:249` — gone with PR 0 anyway) →
  `std::filesystem` via `to_path`; `StringUtils::trim`/`toStr(line)` → `detail::trim` /
  `std::to_string`; `ListUtils::create<std::string>` / `concatenate` (whitelist handling) →
  `detail::split` / `detail::join`. Delete the now-unused includes (`:11–:26`), the vestigial
  `boost/math/fpclassify.hpp` (`:28` — code already uses `std::isnan`), and
  `StringListUtils.h` (`:14` — verify unused, then delete).
* Delete the seeding call `UniqueIdGenerator::setSeed(2453440375)` from `mainInit` (`:77`) — it
  moves to 1.5.
* **Exception translator registry**: `using ExceptionTranslator = bool (*)(std::ostream&);`
  (a translator rethrows `std::current_exception()`, and either prints its recognized type and
  returns `true`, or returns `false`); `registerExceptionTranslator(...)` (small fixed array,
  idempotent); `std::string describeCaughtException()` for the macros. Rework
  `printLastException` (`:682`, called from `END_TEST` `:488` and `:558`): walk translators
  first, then the existing `std::exception` and `...` clauses; delete its
  `BaseException` clause. After this, the `.cpp` names no OpenMS type.

### 1.4 `FuzzyStringComparator.h/.cpp`

* `.h`: drop `Types.h`/`TypeAliases.h`; respell aliases with their definitions — `StringList` →
  `std::vector<std::string>`, `Int` → `int`, `UInt` → `unsigned`. Pure aliases ⇒ source- and
  ABI-identical; `src/topp/FuzzyDiff.cpp` compiles unchanged (verify).
* `.cpp` (`:9–:14`): use the **donor's** version of this file as the guide — it already removed
  `TextFile.h` and most of `SYSTEM/File.h`. On top of that: `File::absolutePath` (six
  report-formatting sites, `:196–:210`) → `std::filesystem::absolute(to_path(...))`;
  file-local copies (unnamed namespace, not installed) of `to_path` (**copy** from
  `SYSTEM/PathUtils.h`, `OPENMS_WINDOWSPLATFORM` → `_WIN32`) and `extractDouble` with its
  helpers `tryParseNaN`/`parseFloat` (**copy** from `StringUtils.cpp:35`/`:229`/definition of
  `StringUtilsHelper::extractDouble`; gate the float-`from_chars` fallback on
  `__cpp_lib_to_chars`, not `OPENMS_NO_FLOAT_FROM_CHARS`).

### 1.5 Test-project support TU (seeding + translator)

New `src/tests/class_tests/openms/source/OpenMSTestSupport.cpp` (not in `executables.cmake`):
a file-local `struct` whose constructor (static init — strictly earlier than today's
`START_TEST`-time seeding) calls `OpenMS::UniqueIdGenerator::setSeed(2453440375)` and
`registerExceptionTranslator(&describeOpenMSException)`, where the translator rethrows, catches
`Exception::BaseException&`, and prints name/file/line/function/`what()` exactly as
`printLastException` does today.

CMake: `openms/CMakeLists.txt:41` → `add_executable(${_class_test} source/${_class_test}.cpp
source/OpenMSTestSupport.cpp)`; `openms_gui/CMakeLists.txt:40` and `:72` add
`${PROJECT_SOURCE_DIR}/../openms/source/OpenMSTestSupport.cpp`. OpenSwathAlgo tests: nothing —
that is the point.

### 1.6 Include-fallout loop

Dropping `StringUtils.h`/`DataValue.h`/`ListUtilsIO.h`/`Types.h` from `ClassTest.h` removes
transitive includes ~720 test TUs may rely on. Measured floor: 26 openms tests need
`DataValue.h`/`ParamValue.h` (23+3) plus 1 `openms_gui` test; tests using
`StringList`/vector-printing may need `TypeAliases.h`/`ListUtilsIO.h` (108 mention the aliases).
Protocol: build, parse missing-declaration errors, add the canonical include at the top of the
failing test, repeat. Budget 30–150 one-line fixes. Never fix by re-adding includes to
`ClassTest.h`.

### 1.7 PR 1 verification

Acceptance 1, 4, 5; plus: `git grep -nE "OpenMS/(DATASTRUCTURES|FORMAT|SYSTEM|CONCEPT/(Types|Exception|PrecisionWrapper|UniqueIdGenerator))" -- src/openms/source/CONCEPT/ClassTest.cpp src/openms/source/CONCEPT/FuzzyStringComparator.cpp src/openms/include/OpenMS/CONCEPT/ClassTest*.h src/openms/include/OpenMS/CONCEPT/FuzzyStringComparator.h`
returns nothing; a deliberately failing `TEST_EQUAL` on a `StringList` and a `DataValue` prints
readable values; a wrong-typed `TEST_EXCEPTION` prints the OpenMS exception name; a
`UniqueId`-bearing format test (featureXML) is byte-stable on Windows and Linux CI.

## PR 2 — the move (pure mechanics)

1. `git mv` the six files to `src/testframework/` (`include/OpenMS/CONCEPT/`: `ClassTest.h`,
   `FuzzyStringComparator.h`, `MacrosTest.h`, `ClassTestUtils.h`; `source/CONCEPT/`: the two
   `.cpp`). Remove their entries from `src/openms/source/CONCEPT/sources.cmake` (`:6`, `:11`)
   and `src/openms/include/OpenMS/CONCEPT/sources.cmake` (`:6`, `:12`, `:20`, + `ClassTestUtils.h`).
2. Strip all `OPENMS_DLLAPI` from the moved headers (50 in `ClassTest.h`, 5 in
   `FuzzyStringComparator.h`, counted on develop) — a static library needs none.
3. `src/testframework/CMakeLists.txt`: take the **donor's** file and delete its
   OpenMS-propagation block (the two `$<TARGET_PROPERTY:OpenMS,...>` stanzas and
   `add_dependencies(... OpenMS)`) — the std-only framework needs no include dirs, definitions,
   or dependency from OpenMS at all. Keep: STATIC, `cxx_std_23`, compiler flags, install of the
   archive under its own component, `install_headers` (paths unchanged ⇒ no test edits),
   `openms_doc_path`. Add `add_subdirectory(testframework)` to `src/CMakeLists.txt` (after
   `openms`; order is irrelevant — no dependency).
4. Link lines:
   * `src/tests/class_tests/openms/CMakeLists.txt:42` →
     `target_link_libraries(${_class_test} OpenMSTestFramework ${OpenMS_LIBRARIES})`;
   * `openms_gui/CMakeLists.txt:41`/`:76` → prepend `OpenMSTestFramework` likewise;
   * `openswathalgo/CMakeLists.txt:39` →
     `target_link_libraries(${i} OpenSwathAlgo OpenMSTestFramework)` — **`OpenMS` is deleted;
     this is the acid test**;
   * `src/topp/CMakeLists.txt`: `target_link_libraries(FuzzyDiff OpenMSTestFramework)` after the
     generic TOPP setup (FuzzyDiff keeps its OpenMS link via TOPPBase).
5. Why no `--as-needed` machinery is needed (put this in a comment where #9919 put its
   `needed`-feature comment): the framework archive references **no** libOpenMS symbols, and
   every openms/openms_gui test binary references libOpenMS through `OpenMSTestSupport.cpp`
   (`setSeed`, `BaseException`) and its own test code — so libOpenMS is a genuine dependency
   wherever it is listed, and link order is unconstrained. The Arrow/Thrift ordering from
   develop (`libOpenMS.so` before the static Arrow archives) is undisturbed because the
   framework introduces no link edge that could reorder it.
6. ABI note for the PR text: libOpenMS stops exporting the ~55 class-test symbols (test-only;
   #9919 measured the set).

Verification: acceptance 2–5; full CI matrix incl. Windows (DLL-free framework archive) and the
wheel builds (`OpenMSTestFramework` stays out of the exported package, as in the donor).

## PR 3 — library-side tidy (any time after PR 1)

1. `Types.h`: delete the `writtenDigits<DataValue>` specialization and its `class DataValue;`
   forward declaration — after PR 1 the only caller was `TEST_REAL_SIMILAR`, which now uses the
   framework's `writtenDigits` (measured: no library code calls it with a `DataValue`).
2. Optional: delete `Internal::OpenMS_locale` (`Types.h:277` decl, `Types.cpp` definition — the
   initializer restores what it sets; net no-op). Fix the two readers:
   `src/topp/OpenSwathFeatureXMLToTSV.cpp:418–:420` (use the literal `"C"`) and the relic
   section in `Types_test.cpp` (delete). `Types.cpp` becomes empty → delete it and its
   `sources.cmake` entry.

## Relationship to #9919

This program subsumes #9919's framework goals from `develop` directly, and by sequencing
std-only *before* extraction it never needs four pieces of #9919's machinery: the user-defined
`needed` link feature, the `LINK_LIBRARY_OVERRIDE_OpenMS` settings, the CMake 3.21→3.24 bump
(which existed only for that feature), and the OpenMS→framework usage-requirement propagation.
`openms_hide_static_archive_symbols()` (#9919's Arrow defense-in-depth) is independent and
optional here, since no link edge is disturbed. If #9919 merges first anyway: PR 0 and the
move-half of PR 2 are already done; apply PR 1 against the moved paths (the previous,
#9919-based revision of this spec — in git history — has those file:line references) and reduce
PR 2 to the link-line simplification (drop the feature uses and the propagation block).

## Risk register

| risk | mitigation |
|---|---|
| numeric-formatter copy drifts from `StringUtils::toStr` | copy `appendNumeric` verbatim; cross-reference comments at both sites; TOPP/FuzzyDiff suite is the oracle |
| `printValue` degrades output where `ListUtilsIO.h` used to be transitive | range branch covers vectors natively; check 1.7's deliberate-failure probes |
| seeding TU dropped from a future test target | `UniqueId`-bearing reference tests fail loudly; comment at all three `add_executable` sites |
| `NDEBUG` vs `OPENMS_ASSERTIONS` divergence in exotic configs | 14 affected files; run one Debug CI job and compare `TEST_PRECONDITION_VIOLATED` outcomes |
| translator unregistered in an exotic binary | `describeCaughtException()` falls back to `std::exception::what()` — degraded output, never wrong results |
| PR 0's opt-in validation misses a *new* test's files | same accepted trade-off as #9919; `VALIDATE_TMP_FILES` documented in the macro docs |

## Explicitly out of scope

`OpenMSCore` / Parts 1–2 of the companion doc (proceeds independently with the libOpenMS split);
folding context into `BaseException::what()` (29 `TEST_EXCEPTION_WITH_MESSAGE` files compare it
verbatim); migrating tests to Catch2 proper.
