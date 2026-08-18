# Implementation spec: std-only `OpenMSTestFramework`

Handoff document for an implementer (human or agent). Self-contained: everything referenced here
was verified on the #9919 branch (`claude/openms-test-framework-refactor-ssutqy`); file:line
references are to that branch. Background and rationale: Part 3 of
[`OpenMSCore-extraction.md`](OpenMSCore-extraction.md).

**Objective.** `OpenMSTestFramework` (`src/testframework/`) compiles and links against the C++
standard library only: zero `#include <OpenMS/...>` outside its own headers, zero undefined
`OpenMS::` symbols in `libOpenMSTestFramework.a` other than its own. OpenMS-specific behavior
(exception naming, ID seeding) is injected from the test projects via registration.

**Acceptance (all four, in order):**
1. `nm -uC OpenMS-build/**/libOpenMSTestFramework.a` lists **no** `OpenMS::` symbols except
   `OpenMS::FuzzyStringComparator::*` / `OpenMS::Internal::ClassTest::*` (its own, defined
   elsewhere in the archive — in practice: no *undefined* `OpenMS::` symbols at all).
2. The five OpenSwathAlgo tests link with `target_link_libraries(${i} OpenSwathAlgo
   OpenMSTestFramework)` — no `OpenMS` anywhere — and pass
   (`readelf -d Datastructures_test` shows no `libOpenMS.so`).
3. Full `ctest` green, including all `TOPP_FuzzyDiff`-based TOPP tests (≈1600 tests exercise the
   fuzzy comparator — they are the oracle for the copied parsing/formatting code).
4. Class-test failure output remains human-equivalent (values still print; OpenMS exception
   names still appear in openms-test failure reports via the translator).

**Ground rules — do NOT:**
* change `BaseException::what()` or `getMessage()` (29 test files compare `what()` verbatim via
  `TEST_EXCEPTION_WITH_MESSAGE`);
* change any TOPP tool, reference file, or library formatting (`StringUtils`, `Types.h` except
  the one deletion in PR B);
* touch anything from Parts 1–2 of the companion doc (`OpenMSCore`, export macros, packaging);
* rename or move the framework's installed headers (install paths are load-bearing for all ~720
  tests).

Work is two PRs. PR A is the substantive one; land it as staged commits in the order given.

---

## PR A — framework goes std-only

Base: the #9919 branch (or `develop` once #9919 merges — `src/testframework/` exists only there).

### A1. New installed header: `src/testframework/include/OpenMS/CONCEPT/ClassTestUtils.h`

Framework-owned utilities needed by `ClassTest.h`'s templates. Namespace
`OpenMS::Internal::ClassTest` (so existing qualified calls keep working); implementation helpers
in a nested `detail` namespace. Contents:

1. **`writtenDigits<T>()`** — the framework's own, replacing the use of `OpenMS::writtenDigits`
   (`Types.h`): floating-point `T` → `std::numeric_limits<T>::digits10`; integral `T` →
   `digits10`; class type convertible to `double` → `std::numeric_limits<double>::digits10`;
   default → 6 (matches the `Types.h` primary template, `Types.h:264`). Needed because the
   OpenSwathAlgo tests use `TEST_REAL_SIMILAR` (3 files) and must not reach `Types.h`.
2. **`isRealType<T>` predicate** — replaces the five overloads at `ClassTest.h:73–113` with one
   function: `template <class T> constexpr bool isRealType(const T&)` returning
   `std::is_floating_point_v<C> || (std::is_class_v<C> && std::is_convertible_v<const C&, double>)`
   with `C = std::remove_cvref_t<T>`. Verified: `DataValue::operator double()`
   (`DataValue.h:175`) and `ParamValue::operator double()` (`ParamValue.h:143`) are implicit, so
   both classify `true`; `std::string`, containers, pointers classify `false`. The
   `TEST_REAL_SIMILAR` macro body (`ClassTest.h:685`) is unchanged — it already calls
   `TEST::isRealType(a)`.
3. **Numeric formatter** `detail::toString(T value)` for arithmetic `T` — **copy** the
   `appendNumeric` implementation from `StringUtils.cpp` (≈`:415`–`:500`: NaN → `"NaN"`, ±inf →
   `"inf"`/`"-inf"`, `std::to_chars` shortest-round-trip with the `1e4`/`1e-2`
   scientific-notation thresholds and the keep-one-fractional-digit rule) plus the plain
   integer `to_chars` path. This must byte-match `StringUtils::toStr` for numeric arguments —
   the string-vs-numeric `TEST_EQUAL` quirk (see A2.3) depends on it, and the ≈720-test suite is
   the oracle. Do not "improve" it.
4. **`detail::printValue(std::ostream&, const T&)`** — printing customization layer:
   if `T` is ostream-streamable (detection idiom) → `os << v`; else if `T` is a range whose
   elements are streamable → print as `[e1, e2, ...]` (element separator `", "`, matching
   `ListUtilsIO.h:29`); else → `os << "<unprintable>"`. Rationale: today vectors print only
   because `ClassTest.h` includes `ListUtilsIO.h`; a framework-local `operator<<` would be
   ambiguous with `OpenMS::operator<<` via ADL for `std::vector<OpenMS-type>`, so a named
   function layer is required, not an operator.
5. **`detail::trim(std::string&)`**, **`detail::split(const std::string&, char)`**,
   **`detail::join(range, const char*)`** — plain implementations (no SIMD). `split` must
   reproduce `ListUtils::create<std::string>` semantics for the `-w` whitelist path (verify
   against `ListUtils.h`'s 13-line implementation when copying; the FuzzyDiff whitelist TOPP
   tests are the oracle).

### A2. `ClassTest.h` edits

1. Delete includes (`:36–:49`): `PrecisionWrapper.h`, `Types.h`, `StringUtils.h`, `DataValue.h`,
   `OpenMSConfig.h`, `config.h`, `ListUtilsIO.h`. Add `ClassTestUtils.h`. (`MacrosTest.h` stays.)
2. Delete the five `isRealType` overloads (`:73–:113`) — replaced by A1.2.
3. `testEqual` (`:288`): replace `StringUtils::toStr(expression_2)` (`:303`) with
   `detail::toString(expression_2)`; replace the three bare `stdcout << expression_N` prints
   (`:314–:329`, and the analogous `testNotEqual` at `:402–:408`) with
   `detail::printValue(stdcout, expression_N)`. Keep the existing `XMLCh*` and enum
   `if constexpr` branches exactly as they are.
4. `TEST_REAL_SIMILAR` (`:685`): change `writtenDigits(a)` → `TEST::writtenDigits(a)` (twice) so
   it binds to A1.1 instead of requiring `Types.h` in the test TU.
5. Exception macros — make the macro text std-only. In `TEST_EXCEPTION` (and `TEST_EXCEPTION_WITH_MESSAGE`,
   and any sibling like `TEST_PRECONDITION_VIOLATED`'s expansion path):
   * keep `catch (exception_type&)` (and the `what()` comparison clause in `_WITH_MESSAGE` —
     unchanged, it compares `et.what()` verbatim);
   * **delete** the `catch (::OpenMS::Exception::BaseException& e)` clause
     (`TEST_EXCEPTION` body, and `_WITH_MESSAGE` at `:943`ff);
   * collapse the remaining fallbacks to `catch (...) { TEST::exception = 2;
     TEST::exception_name = TEST::describeCaughtException(); }` where
     `describeCaughtException()` is a new `.cpp` function (A3.2). Adjust the switch that prints
     the outcome accordingly (the "wrong exception thrown" message now prints
     `exception_name` from the describer — richer than today for openms tests, `what()`-based
     for everything else).
6. Precondition macros (`:903`/`:920`): `#ifdef OPENMS_ASSERTIONS` → `#ifndef NDEBUG`.
   Semantics note: today `OPENMS_ASSERTIONS` is set from `CMAKE_BUILD_TYPE STREQUAL Debug`
   (never true for multi-config Debug); `!NDEBUG` tracks the actual compile mode of the test TU,
   which is what must agree with how the *library*'s `OPENMS_PRECONDITION` was compiled.
   14 test files use `TEST_PRECONDITION_VIOLATED` — confirm they still pass in a Debug build.
7. Replace any remaining `OPENMS_WINDOWSPLATFORM` in framework files with `_WIN32`
   (one use at `ClassTest.h:76` region; also check `.cpp`s).

### A3. `ClassTest.cpp` edits

1. Delete includes: `UniqueIdGenerator.h`, `ListUtils.h` (and transitively everything OpenMS).
   Replace: `UniqueIdGenerator::setSeed(2453440375)` (`:78`) → **delete** (moves to A5);
   `ListUtils::create<std::string>(whitelist_, ',')` (`:210`) → `detail::split`;
   `ListUtils::concatenate(...)` (`:216`) → `detail::join`; `StringUtils::trim` (`:113–:114`) →
   `detail::trim`; `StringUtils::toStr(line)` (`:250`) → `std::to_string(line)`.
2. **Exception translator registry.** Add to the framework:
   ```cpp
   using ExceptionTranslator = bool (*)(std::ostream&); // rethrows current_exception; prints and
                                                        // returns true if it recognized the type
   void registerExceptionTranslator(ExceptionTranslator);   // idempotent, small fixed array
   std::string describeCaughtException();                   // for the macros (A2.5)
   ```
   Rework `printLastException` (currently rethrows and catches `BaseException` explicitly,
   printing name/file/line/function/what): first walk registered translators (each rethrows and
   returns false if the type is not theirs), then the existing `std::exception` and `...`
   clauses. `describeCaughtException()` is the string-returning analogue for the macro path.
   Net effect: the `.cpp` no longer mentions any OpenMS type; openms tests keep today's
   diagnostics via A5.
3. `FuzzyStringComparator.cpp`: replace `#include <OpenMS/DATASTRUCTURES/StringUtils.h>` and
   `#include <OpenMS/SYSTEM/PathUtils.h>` with file-local copies in an unnamed namespace:
   `extractDouble` + its helpers `tryParseNaN`/`parseFloat` (**copy** from
   `StringUtils.cpp:382` region; gate the float-`from_chars` fallback on
   `__cpp_lib_to_chars`, not `OPENMS_NO_FLOAT_FROM_CHARS`), and `to_path` (**copy** the 20-line
   header-only function from `SYSTEM/PathUtils.h`, `OPENMS_WINDOWSPLATFORM` → `_WIN32`).
4. `FuzzyStringComparator.h`: drop `Types.h`/`TypeAliases.h` includes; respell the aliases with
   their underlying types — `StringList` → `std::vector<std::string>`, `Int` → `int`, `UInt` →
   `unsigned` (`:93–:99` `get/setWhitelist`, `:341–:343` members, verbose-level setters). These
   are pure aliases, so the API is source- and ABI-identical; `src/topp/FuzzyDiff.cpp`
   (`:110–:132`) needs no change — verify by compiling it.

### A4. Framework CMake decoupling — `src/testframework/CMakeLists.txt`

Delete the two `target_include_directories(... $<TARGET_PROPERTY:OpenMS,...>)` blocks, the
`target_compile_definitions(... $<TARGET_PROPERTY:OpenMS,...>)`, and
`add_dependencies(OpenMSTestFramework OpenMS)` (`:49–:61`). Add `ClassTestUtils.h` to
`OpenMSTestFramework_sources_h`. Rewrite the header comment: the framework depends on the C++
standard library only. Everything else (static, install, components) stays.

### A5. Test-project support TU (seeding + translator)

New file `src/tests/class_tests/openms/source/OpenMSTestSupport.cpp` (not a test; not listed in
`executables.cmake`):

```cpp
// Registers OpenMS-specific behavior with the std-only class-test framework:
// deterministic unique IDs for reference-file stability, and rich reporting
// for unexpected OpenMS exceptions.
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
namespace {
bool describeOpenMSException(std::ostream& os) { /* rethrow; catch BaseException&;
  print name/file/line/function/what exactly as printLastException does today; return true */ }
struct OpenMSTestSupport {
  OpenMSTestSupport() {
    OpenMS::UniqueIdGenerator::setSeed(2453440375);
    OpenMS::Internal::ClassTest::registerExceptionTranslator(&describeOpenMSException);
  }
} openms_test_support_init;
}
```

Seeding at static init runs strictly before `START_TEST` did — same or better determinism.

CMake: `src/tests/class_tests/openms/CMakeLists.txt:49` →
`add_executable(${_class_test} source/${_class_test}.cpp source/OpenMSTestSupport.cpp)`;
`openms_gui/CMakeLists.txt:46` and `:89` reference the same file via
`${PROJECT_SOURCE_DIR}/../openms/source/OpenMSTestSupport.cpp`. OpenSwathAlgo tests: nothing —
that is the point.

### A6. OpenSwathAlgo test link lines

`src/tests/class_tests/openswathalgo/CMakeLists.txt`: `target_link_libraries(${i} OpenSwathAlgo
OpenMSTestFramework)`; delete `"$<LINK_LIBRARY:needed,OpenMS>"` and its four-line comment. The
openms/openms_gui test projects keep `"$<LINK_LIBRARY:needed,OpenMS>"` and
`LINK_LIBRARY_OVERRIDE_OpenMS` unchanged (libOpenMS is their subject); update the comment at
`openms/CMakeLists.txt:55–:61` to say the framework itself no longer needs OpenMS.

### A7. Include-fallout loop (the long tail)

Dropping `StringUtils.h`/`DataValue.h`/`ListUtilsIO.h`/`Types.h` from `ClassTest.h` removes
transitive includes ~720 test TUs may rely on. Known floor: 26 openms tests need
`DataValue.h`/`ParamValue.h` (23 + 3, measured) and 1 `openms_gui` test needs `DataValue.h`;
tests using `StringList` (108 mention it) or printing vectors may need
`DATASTRUCTURES/ListUtilsIO.h` / `TypeAliases.h`. Protocol: build `ctest` targets, parse
missing-declaration errors, insert the canonical include at the top of the failing test file,
repeat until clean. Every fix is one line; budget 30–150 files; do not "fix" by re-adding
includes to `ClassTest.h`.

### A8. Verification for PR A

1. Acceptance checks 1–4 above (in that order — the `nm` check first, it is cheap and precise).
2. `TOPP_FuzzyDiff` suite green (parsing copy), plus one manual spot check: run `FuzzyDiff` on a
   pair of files with numbers near the `1e4`/`1e-2` thresholds and NaN/inf tokens; identical
   verdicts before/after.
3. A deliberately failing `TEST_EQUAL` on a `StringList` and on a `DataValue` prints readable
   values (printValue path), and a deliberately wrong-typed `TEST_EXCEPTION` in an openms test
   prints the OpenMS exception *name* (translator path).
4. Windows + macOS CI: the seeding TU must show identical reference-file behavior (any
   `UniqueId`-bearing format test, e.g. featureXML round-trips).
5. `git grep -n "OpenMS/" src/testframework/` returns only the framework's own headers.

CHANGELOG entry; update `AGENTS.md` if it describes the framework's dependencies.

---

## PR B — library-side tidy (after A)

1. `Types.h`: delete the `writtenDigits<DataValue>` specialization and the `class DataValue;`
   forward declaration (`:222–:229`). Verified unreferenced by library code once
   `TEST_REAL_SIMILAR` calls the framework's `writtenDigits` (measured: no library caller
   passes a `DataValue` today).
2. Optional, same PR: delete `Internal::OpenMS_locale` (`Types.cpp`, declaration at
   `Types.h:277`) — the initializer restores what it sets (net no-op) and the global holds
   `"C"`. Fix the two readers: `src/topp/OpenSwathFeatureXMLToTSV.cpp:418–:420` (use the
   literal) and the relic section in `Types_test.cpp:27–:39` (delete). `Types.cpp` becomes
   empty → delete the file and its `sources.cmake` entry.

## Explicitly out of scope

* `OpenMSCore` / everything in Parts 1–2 of the companion doc (proceeds independently, driven by
  the libOpenMS split).
* Folding file/line/name into `BaseException::what()` — blocked by the 29
  `TEST_EXCEPTION_WITH_MESSAGE` files comparing `what()` verbatim; revisit separately if wanted.
* Migrating tests to Catch2 proper.

## Risk register

| risk | mitigation |
|---|---|
| numeric-formatter copy drifts from `StringUtils::toStr` | copy `appendNumeric` verbatim, do not modify; suite is the oracle; keep a comment naming the source of truth in both places |
| `printValue` falls back to `<unprintable>` where today's code compiled via `ListUtilsIO.h` | the range branch covers vectors natively; anything else that regresses shows up in check A8.3 |
| seeding TU dropped from a future test target | the featureXML-class reference tests fail loudly on any unseeded run; comment at the `add_executable` sites |
| `NDEBUG` vs `OPENMS_ASSERTIONS` divergence in exotic configs | 14 affected files; run one Debug-config CI job and compare `TEST_PRECONDITION_VIOLATED` outcomes |
| translator not registered in some exotic test binary | `describeCaughtException()` falls back to `std::exception::what()` — degraded output, never wrong results |
