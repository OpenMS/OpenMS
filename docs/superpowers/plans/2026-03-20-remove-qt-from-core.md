# Remove Qt from OpenMS Core Library — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove Qt6::Core as a dependency of `libOpenMS` so pyOpenMS and other consumers can build without Qt.

**Architecture:** Remove `String(const QString&)` constructor and `toQString()` from the core `String` class. Provide `toQString()`/`fromQString()` free functions in `Qt5Port.h` (in `openms_gui`). Migrate ~60 GUI files to free functions, ~12 TOPP tools to inline conversions or free functions, remove ~17 stale Qt includes, and add explicit `Qt6::Core` link deps where needed. Update CMake config so Qt is only required for GUI/TOPP builds, not core lib consumers.

**Tech Stack:** C++ (std::string, QString), CMake, Qt6::Core

**Spec:** `docs/superpowers/specs/2026-03-20-remove-qt-from-core-lib-design.md`

---

### Task 1: Add compatibility free functions to Qt5Port.h

**Files:**
- Modify: `src/openms_gui/include/OpenMS/VISUAL/MISC/Qt5Port.h`

This must land first — all subsequent migrations depend on these functions being available.

- [ ] **Step 1: Add `#include <QString>` and free functions to Qt5Port.h**

Add before the closing `} // namespace OpenMS`:

```cpp
#include <QString>

/// Convert OpenMS::String to QString (replaces String::toQString())
inline QString toQString(const String& s) { return QString::fromStdString(s); }

/// Construct OpenMS::String from QString (replaces String(const QString&))
inline String fromQString(const QString& s) { return String(s.toStdString()); }
```

- [ ] **Step 2: Build to verify no compilation errors**

Run: `cmake --build OpenMS-build --target OpenMS_GUI -j$(nproc)`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/openms_gui/include/OpenMS/VISUAL/MISC/Qt5Port.h
git commit -m "feat: add toQString/fromQString free functions to Qt5Port.h"
```

---

### Task 2: Migrate GUI library files to free functions

**Files:**
- Modify: ~60 files in `src/openms_gui/source/` and ~7 headers in `src/openms_gui/include/`

Mechanical search-and-replace across the GUI lib. Every `.toQString()` member call becomes a `toQString()` free function call. Every `String(some_qstring)` implicit conversion becomes `fromQString(some_qstring)`. Add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>` where not already present.

- [ ] **Step 1: Find all `.toQString()` calls in `src/openms_gui/`**

Run: `grep -rn '\.toQString()' src/openms_gui/ --include='*.cpp' --include='*.h' | wc -l`

Record the count for verification later.

- [ ] **Step 2: Replace `.toQString()` calls with free function**

For each file, replace patterns:
- `expr.toQString()` → `toQString(expr)` (for simple expressions)
- `String(...).toQString()` → `QString::fromStdString(String(...))` or `toQString(String(...))`

Add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>` at the top of each modified file if not already present.

Key headers to check (have inline usage):
- `src/openms_gui/include/OpenMS/VISUAL/FileWatcher.h` (lines 62, 68)

- [ ] **Step 3: Replace `String(qstring_expr)` conversions with `fromQString()`**

Search for patterns like `String(some_expr)` where `some_expr` returns a `QString`. These are harder to grep for mechanically — look for:
- `String(QDir(...)...)`
- `String(QFileInfo(...)...)`
- `String(someQString)`
- Assignments: `String s = qstring_expr;`

Replace with `fromQString(qstring_expr)` or `qstring_expr.toStdString()`.

- [ ] **Step 4: Build to verify compilation**

Run: `cmake --build OpenMS-build --target OpenMS_GUI -j$(nproc)`
Expected: PASS

- [ ] **Step 5: Verify no remaining `.toQString()` member calls**

Run: `grep -rn '\.toQString()' src/openms_gui/ --include='*.cpp' --include='*.h'`
Expected: Only `Qt5Port.h` itself (the free function definition) should remain.

- [ ] **Step 6: Commit**

```bash
git add src/openms_gui/
git commit -m "refactor: migrate GUI lib from String::toQString() to free functions"
```

---

### Task 3: Migrate TOPP tools with GUI link to free functions

**Files:**
- Modify: `src/topp/ExecutePipeline.cpp`, `src/topp/ImageCreator.cpp`, `src/topp/INIUpdater.cpp`

These 3 tools link `OpenMS_GUI` (via `TOPP_executables_with_GUIlib` in `src/topp/executables.cmake:163-170`) so they can use `Qt5Port.h`.

- [ ] **Step 1: Migrate ExecutePipeline.cpp**

Add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>` and replace:
- `str.toQString()` → `toQString(str)`
- `String(qstr)` → `fromQString(qstr)`

- [ ] **Step 2: Migrate ImageCreator.cpp**

Same pattern.

- [ ] **Step 3: Migrate INIUpdater.cpp**

Same pattern.

- [ ] **Step 4: Build to verify**

Run: `cmake --build OpenMS-build --target ExecutePipeline ImageCreator INIUpdater -j$(nproc)`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/topp/ExecutePipeline.cpp src/topp/ImageCreator.cpp src/topp/INIUpdater.cpp
git commit -m "refactor: migrate GUI-linked TOPP tools to Qt5Port.h free functions"
```

---

### Task 4: Migrate TOPP tools without GUI link to inline conversions

**Files:**
- Modify: `src/topp/MetaProSIP.cpp`, `src/topp/GenericWrapper.cpp`, `src/topp/IDRipper.cpp`, `src/topp/MzMLSplitter.cpp`, `src/topp/CometAdapter.cpp`, `src/topp/OpenSwathFileSplitter.cpp`, `src/topp/AssayGeneratorMetaboSirius.cpp`, `src/topp/OpenSwathWorkflow.cpp`, `src/topp/MaRaClusterAdapter.cpp`

These tools link only `OpenMS` (not `OpenMS_GUI`), so they cannot include `Qt5Port.h`. Inline the conversions directly.

- [ ] **Step 1: For each file, replace `.toQString()` with `QString::fromStdString()`**

Pattern: `str.toQString()` → `QString::fromStdString(str)`

**Watch out for round-trip conversions:** If `str.toQString()` is passed to a function expecting `String` or `const String&`, the `.toQString()` is a no-op round-trip — just remove it entirely and pass `str` directly. Example in `MaRaClusterAdapter.cpp:369`:
```cpp
// Before (round-trip: String→QString→String via implicit constructor):
bool copy_status = File::copyDirRecursively(tmp_dir.getPath().toQString(), maracluster_output_directory.toQString());
// After (just remove the .toQString() calls):
bool copy_status = File::copyDirRecursively(tmp_dir.getPath(), maracluster_output_directory);
```

- [ ] **Step 2: Replace `String(qstring_expr)` with `.toStdString()`**

Pattern: `String(qdir.absolutePath())` → `qdir.absolutePath().toStdString()`

Key site: `src/topp/OpenSwathWorkflow.cpp:668`:
```cpp
// Before:
String tmp_dir = String(QDir(getStringOption_("tempDirectory").c_str()).absolutePath()).ensureLastChar('/');
// After:
String tmp_dir = String(QDir(getStringOption_("tempDirectory").c_str()).absolutePath().toStdString()).ensureLastChar('/');
```

- [ ] **Step 3: Build to verify**

Run: `cmake --build OpenMS-build --target MetaProSIP GenericWrapper IDRipper MzMLSplitter CometAdapter OpenSwathFileSplitter AssayGeneratorMetaboSirius OpenSwathWorkflow MaRaClusterAdapter -j$(nproc)`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add src/topp/MetaProSIP.cpp src/topp/GenericWrapper.cpp src/topp/IDRipper.cpp src/topp/MzMLSplitter.cpp src/topp/CometAdapter.cpp src/topp/OpenSwathFileSplitter.cpp src/topp/AssayGeneratorMetaboSirius.cpp src/topp/OpenSwathWorkflow.cpp src/topp/MaRaClusterAdapter.cpp
git commit -m "refactor: migrate non-GUI TOPP tools to inline Qt string conversions"
```

---

### Task 5: Remove stale Qt includes from TOPP tools

**Files:**
- Modify: `src/topp/PSMFeatureExtractor.cpp`, `src/topp/NucleicAcidSearchEngine.cpp`, `src/topp/MapStatistics.cpp`, `src/topp/FeatureFinderMultiplex.cpp`, `src/topp/LuciphorAdapter.cpp`, `src/topp/QCExporter.cpp`, `src/topp/QCImporter.cpp`, `src/topp/QCMerger.cpp`, `src/topp/Resampler.cpp`

Also remove the stale `<QtCore/QProcess>` include from `MaRaClusterAdapter.cpp` (line 28 — its `.toQString()` calls were handled in Task 4).

- [ ] **Step 1: Remove unused Qt includes from each file**

Remove lines like:
- `#include <QtCore/QFile>`
- `#include <QtCore/QDir>`
- `#include <QtCore/QProcess>`
- `#include <QtCore/QProcessEnvironment>`
- `#include <QtCore/QString>`
- `#include <QtGui/QImage>` (Resampler.cpp — note: this is a QtGui include, not QtCore)

Verify each is truly unused by checking no Qt symbols remain in the file after the include removal.

- [ ] **Step 2: Build to verify**

Run: `cmake --build OpenMS-build --target PSMFeatureExtractor NucleicAcidSearchEngine MapStatistics FeatureFinderMultiplex LuciphorAdapter QCExporter QCImporter QCMerger Resampler MaRaClusterAdapter -j$(nproc)`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/topp/PSMFeatureExtractor.cpp src/topp/NucleicAcidSearchEngine.cpp src/topp/MapStatistics.cpp src/topp/FeatureFinderMultiplex.cpp src/topp/LuciphorAdapter.cpp src/topp/QCExporter.cpp src/topp/QCImporter.cpp src/topp/QCMerger.cpp src/topp/Resampler.cpp src/topp/MaRaClusterAdapter.cpp
git commit -m "refactor: remove stale Qt includes from TOPP tools"
```

---

### Task 6: Migrate and clean up test files

**Files:**
- Modify: `src/tests/class_tests/openms/source/String_test.cpp`
- Modify: `src/tests/class_tests/openms/source/MzMLSqliteHandler_test.cpp`
- Modify: `src/tests/class_tests/openms/source/PythonInfo_test.cpp`
- Modify: `src/tests/class_tests/openms/source/ToolDescriptionFile_test.cpp`

Also remove stale includes from:
- `src/tests/class_tests/openms/source/MzMLSqliteSwathHandler_test.cpp`
- `src/tests/class_tests/openms/source/StringListUtils_test.cpp`
- `src/tests/class_tests/openms/source/SqMassFile_test.cpp`
- `src/tests/class_tests/openms/source/OPXLHelper_test.cpp`
- `src/tests/class_tests/openms/source/DataValue_test.cpp`
- `src/tests/class_tests/openms/source/SiriusFragmentAnnotation_test.cpp`
- `src/tests/class_tests/openms/source/XQuestResultXMLFile_test.cpp`

- [ ] **Step 1: Remove Qt-specific test cases from String_test.cpp**

Remove test cases that test `String(const QString&)` constructor and `toQString()` method. Search for `QString` in the test file and remove the relevant `START_SECTION`/`END_SECTION` blocks.

- [ ] **Step 2: Migrate MzMLSqliteHandler_test.cpp, PythonInfo_test.cpp, ToolDescriptionFile_test.cpp**

Replace `.toQString()` calls with `QString::fromStdString(str)` and `String(qstr)` with `qstr.toStdString()`. These tests don't link `openms_gui`.

- [ ] **Step 3: Remove stale Qt includes from 7 test files**

Remove unused `#include <QtCore/Q...>` lines from each file listed above.

- [ ] **Step 4: Build and run affected tests**

Run: `cmake --build OpenMS-build -j$(nproc) && ctest --test-dir OpenMS-build -R "String_test|MzMLSqliteHandler_test|PythonInfo_test|ToolDescriptionFile_test" -V`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/tests/class_tests/openms/source/
git commit -m "refactor: remove Qt dependencies from core library tests"
```

---

### Task 7: Remove Qt from core library (String.h, String.cpp, TOPPBase.cpp, CMakeLists.txt)

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/String.h`
- Modify: `src/openms/source/DATASTRUCTURES/String.cpp`
- Modify: `src/openms/source/APPLICATIONS/TOPPBase.cpp`
- Modify: `src/openms/CMakeLists.txt`

This is the core change. All callers must be migrated (Tasks 2-6) before this step.

- [ ] **Step 1: Remove Qt from String.h**

Remove:
- Line 18: `class QString;` forward declaration
- Line 73: `String(const QString& s);` constructor declaration
- Line 335: `QString toQString() const;` method declaration

- [ ] **Step 2: Remove Qt from String.cpp**

Remove:
- Line 17: `#include <QString>`
- Lines 40-43: `String(const QString&)` constructor implementation
- Lines 305-308: `toQString()` method implementation

- [ ] **Step 3: Remove stale QDateTime include from TOPPBase.cpp**

Remove line 46: `#include <QtCore/QDateTime>`

- [ ] **Step 4: Remove Qt6::Core from CMakeLists.txt**

In `src/openms/CMakeLists.txt` line 61, remove `Qt6::Core` from the `OPENMS_DEP_LIBRARIES` list.

- [ ] **Step 5: Build core library**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc)`
Expected: PASS — core lib compiles with zero Qt.

- [ ] **Step 6: Commit**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/String.h src/openms/source/DATASTRUCTURES/String.cpp src/openms/source/APPLICATIONS/TOPPBase.cpp src/openms/CMakeLists.txt
git commit -m "refactor: remove Qt dependency from OpenMS core library"
```

---

### Task 8: Add explicit Qt6::Core link deps for TOPP tools and tests

**Files:**
- Modify: `src/topp/CMakeLists.txt`
- Modify: `src/tests/class_tests/openms/CMakeLists.txt`

After removing Qt6::Core from libOpenMS's PUBLIC deps, targets that directly use Qt headers need explicit linking.

- [ ] **Step 1: Add Qt6::Core to affected TOPP tools in `src/topp/CMakeLists.txt`**

Add after line 71 (the existing `target_link_libraries` block):

```cmake
# Qt6::Core for tools that directly use Qt headers
target_link_libraries(AssayGeneratorMetaboSirius Qt6::Core)
target_link_libraries(CometAdapter Qt6::Core)
target_link_libraries(GenericWrapper Qt6::Core)
target_link_libraries(IDRipper Qt6::Core)
target_link_libraries(MetaProSIP Qt6::Core)
target_link_libraries(MzMLSplitter Qt6::Core)
target_link_libraries(OpenSwathFileSplitter Qt6::Core)
target_link_libraries(OpenSwathWorkflow Qt6::Core)
target_link_libraries(QCEmbedder Qt6::Core)
target_link_libraries(QCExtractor Qt6::Core)
target_link_libraries(QCShrinker Qt6::Core)
```

- [ ] **Step 2: Add Qt6::Core to affected tests in `src/tests/class_tests/openms/CMakeLists.txt`**

Add after line 74 (the existing special link deps):

```cmake
# Qt6::Core for tests that directly use Qt headers
target_link_libraries(FIAMSScheduler_test Qt6::Core)
target_link_libraries(MzMLSqliteHandler_test Qt6::Core)
target_link_libraries(PeptideIndexing_test Qt6::Core)
target_link_libraries(PythonInfo_test Qt6::Core)
target_link_libraries(ToolDescriptionFile_test Qt6::Core)
```

- [ ] **Step 3: Full build to verify everything compiles**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: PASS — full build succeeds.

- [ ] **Step 4: Commit**

```bash
git add src/topp/CMakeLists.txt src/tests/class_tests/openms/CMakeLists.txt
git commit -m "build: add explicit Qt6::Core link deps for TOPP tools and tests"
```

---

### Task 9: Update CMake configuration for Qt-free core lib

**Files:**
- Modify: `cmake/cmake_findExternalLibs.cmake`
- Modify: `cmake/OpenMSConfig.cmake.in`
- Modify: `cmake/add_library_macros.cmake`

- [ ] **Step 1: Keep Qt Core find_package but decouple from core lib in `cmake/cmake_findExternalLibs.cmake`**

**Important:** `find_package(Qt6 Core)` must stay **outside** `if(WITH_GUI)` because non-GUI TOPP tools (MetaProSIP, GenericWrapper, etc.) directly use Qt headers and need `Qt6::Core` linked. These tools are built unconditionally (not gated by `WITH_GUI`).

Lines 290-291: Keep the `find_package` call but update the variable name and comment to clarify it's no longer for the core lib:

Before:
```cmake
set(OpenMS_QT_COMPONENTS Core CACHE INTERNAL "QT components for core lib")
find_package(Qt6 ${QT_MIN_VERSION} COMPONENTS ${OpenMS_QT_COMPONENTS} REQUIRED)
```

After:
```cmake
# Qt6::Core is needed by TOPP tools and GUI, but NOT by libOpenMS itself
find_package(Qt6 ${QT_MIN_VERSION} COMPONENTS Core REQUIRED)
```

Remove the `OpenMS_QT_COMPONENTS` cache variable since it's no longer used by the core lib (it was only referenced in `OpenMSConfig.cmake.in`).

- [ ] **Step 2: Update `cmake/OpenMSConfig.cmake.in`**

Line 25: `find_dependency(Qt6 @QT_MIN_VERSION@ COMPONENTS @OpenMS_QT_COMPONENTS@)`

Remove this line. `libOpenMS` no longer needs Qt, so downstream consumers of only `libOpenMS` should not be forced to find Qt. The `@OpenMS_QT_COMPONENTS@` variable has been removed in Step 1.

Note: Consumers using `OpenMS_GUI` get Qt via its own dependencies. If a downstream project uses TOPP tools directly, they need to find Qt themselves.

- [ ] **Step 3: Disable AUTOMOC for core lib in `cmake/add_library_macros.cmake`**

Line 134 sets `AUTOMOC ON` for all libraries. The core lib no longer uses any Qt metaobject features. Change to only enable AUTOMOC when Qt is available:

```cmake
if(TARGET Qt6::moc)
  set_target_properties(${openms_add_library_TARGET_NAME} PROPERTIES AUTOMOC ON)
endif()
```

Or pass a parameter to `openms_add_library` to control AUTOMOC per-target.

- [ ] **Step 4: Reconfigure and full build**

Run: `cmake OpenMS-build && cmake --build OpenMS-build -j$(nproc)`
Expected: PASS (re-runs CMake using existing cache, then builds)

- [ ] **Step 5: Commit**

```bash
git add cmake/cmake_findExternalLibs.cmake cmake/OpenMSConfig.cmake.in cmake/add_library_macros.cmake
git commit -m "build: make Qt optional for core lib in CMake configuration"
```

---

### Task 10: Verify pyOpenMS builds without Qt and run full test suite

**Files:** None (verification only)

- [ ] **Step 1: Build pyOpenMS**

Run: `cmake --build OpenMS-build --target pyopenms -j$(nproc)`
Expected: PASS — pyOpenMS builds without any Qt symbols.

- [ ] **Step 2: Run pyOpenMS tests**

Run: `cd /tmp && PYTHONPATH=/home/sachsenb/Development/tmp2/OpenMS/OpenMS-build/pyOpenMS python3 -m pytest /home/sachsenb/Development/tmp2/OpenMS/src/pyOpenMS/tests/ -v`
Expected: PASS

- [ ] **Step 3: Run full C++ test suite**

Run: `ctest --test-dir OpenMS-build -j$(nproc)`
Expected: All tests pass.

- [ ] **Step 4: Verify no Qt symbols in libOpenMS**

Run: `nm -D OpenMS-build/lib/libOpenMS.so | grep -i qstring | head -20`
Expected: No Qt string symbols found.

- [ ] **Step 5: Commit any remaining fixes if needed**

---

## Task Dependency Graph

```
Task 1 (Qt5Port.h free functions)
  ├── Task 2 (GUI lib migration) ──┐
  ├── Task 3 (TOPP GUI tools)     ├── Task 7 (Remove Qt from core)
  ├── Task 4 (TOPP non-GUI tools) │     ├── Task 8 (Explicit Qt links)
  ├── Task 5 (Stale TOPP includes)│     ├── Task 9 (CMake config)
  └── Task 6 (Tests)──────────────┘     └── Task 10 (Verification)
```

Tasks 2-6 can be parallelized (independent migrations). Task 7 depends on all of 2-6. Tasks 8-9 depend on 7. Task 10 depends on all.
