# Remove Qt Dependency from OpenMS Core Library

**Date:** 2026-03-20
**Goal:** Make `libOpenMS` buildable without Qt so that pyOpenMS (and other consumers) do not require Qt as a transitive dependency.

## Background

Qt6::Core is currently a PUBLIC dependency of `libOpenMS`. After a series of refactors that replaced QDir/QFile (`std::filesystem`), QProcess (`boost::process`), QDateTime (`std::chrono`), and QJson (`nlohmann/json`), the only remaining Qt usage in the core library is String/QString interop in 3 files.

### Remaining Qt Usage (pre-change)

| File | Qt Usage |
|------|----------|
| `src/openms/include/OpenMS/DATASTRUCTURES/String.h` | Forward decl `class QString`, `String(const QString&)` constructor, `toQString()` method |
| `src/openms/source/DATASTRUCTURES/String.cpp` | `#include <QString>`, constructor and method implementations |
| `src/openms/source/APPLICATIONS/TOPPBase.cpp` | `#include <QtCore/QDateTime>` (stale, unused) |
| `src/openms/CMakeLists.txt` | `Qt6::Core` in `OPENMS_DEP_LIBRARIES` (PUBLIC) |

No Q_OBJECT, signals/slots, or Qt metaobject system usage exists in the core library.

## Approach

Minimal surgical removal: remove Qt from the core lib's public API and provide replacements for callers:
- **GUI lib callers** (`openms_gui`): use free functions in the existing `Qt5Port.h` compatibility header
- **Non-GUI callers** (TOPP tools, tests): inline the conversions directly (`QString::fromStdString(str)` / `qstr.toStdString()`) since they cannot include headers from `openms_gui`

## Changes

### 1. Core Library (`src/openms/`)

**`String.h`** — Remove:
- `class QString;` forward declaration (line 18)
- `String(const QString& s);` constructor (line 73)
- `QString toQString() const;` method (line 335)

**`String.cpp`** — Remove:
- `#include <QString>` (line 17)
- `String(const QString&)` constructor implementation (lines 40-43)
- `toQString()` method implementation (lines 305-308)

**`TOPPBase.cpp`** — Remove:
- `#include <QtCore/QDateTime>` (line 46, stale/unused — the source comment says "still needed" but no QDateTime symbols are used in the file)

**`CMakeLists.txt`** — Remove:
- `Qt6::Core` from `OPENMS_DEP_LIBRARIES` (line 61)

### 2. CMake Configuration Changes

Removing `Qt6::Core` from `OPENMS_DEP_LIBRARIES` alone is not sufficient. The build system also references Qt in several other places that must be updated:

**`cmake/cmake_findExternalLibs.cmake`** (line 291-299):
- The `find_package(Qt6 ... COMPONENTS Core REQUIRED)` call and the `OpenMS_QT_COMPONENTS` variable are used to find Qt for the core lib. Make Qt6::Core conditional — only required when building `openms_gui` or TOPP tools, not when building only `libOpenMS`.
- Move the Qt Core find_package into the GUI section (alongside the existing `find_package(Qt6 ... ${OpenMS_GUI_QT_COMPONENTS})` block at line 330).

**`cmake/OpenMSConfig.cmake.in`** (line 25):
- `find_dependency(Qt6 ... COMPONENTS @OpenMS_QT_COMPONENTS@)` forces downstream consumers to have Qt. Remove Qt Core from this list so that consumers of only `libOpenMS` do not need Qt. Qt should only be required if the consumer also uses `OpenMS_GUI`.

**`cmake/add_library_macros.cmake`** (line 134):
- `AUTOMOC ON` is set for all libraries including `libOpenMS`. Since the core lib no longer uses any Qt metaobject features (no Q_OBJECT, no signals/slots), disable AUTOMOC for the `OpenMS` target specifically, or make it conditional on GUI builds only.

### 3. Compatibility Layer (`src/openms_gui/`)

**`src/openms_gui/include/OpenMS/VISUAL/MISC/Qt5Port.h`** — Add two free functions and an explicit `<QString>` include:

```cpp
#include <QString>

/// Convert OpenMS::String to QString
inline QString toQString(const String& s) { return QString::fromStdString(s); }

/// Construct OpenMS::String from QString
inline String fromQString(const QString& s) { return s.toStdString(); }
```

Note: Add `#include <QString>` explicitly rather than relying on `<QStringList>` transitively.

### 4. GUI Library Migration (`src/openms_gui/`)

~60 source files and headers that call `.toQString()` or use `String(qstr)` implicit conversion:
- Replace `str.toQString()` with `toQString(str)` (free function from `Qt5Port.h`)
- Replace `String(some_qstring)` with `fromQString(some_qstring)` (free function from `Qt5Port.h`)
- Add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>` where not already included

This includes both `.cpp` files and headers with inline usage (e.g., `FileWatcher.h` lines 62, 68).

These files already link `OpenMS_GUI` and have the include path available.

### 5. TOPP Tool Migration (`src/topp/`)

12 files that call `.toQString()` or use `String(qstr)`:

**Link openms_gui (can use `Qt5Port.h`):** `ExecutePipeline.cpp`, `ImageCreator.cpp`, `INIUpdater.cpp`

**Link only OpenMS (inline conversions directly):** `MetaProSIP.cpp`, `GenericWrapper.cpp`, `IDRipper.cpp`, `MzMLSplitter.cpp`, `CometAdapter.cpp`, `OpenSwathFileSplitter.cpp`, `AssayGeneratorMetaboSirius.cpp`, `OpenSwathWorkflow.cpp` (line 668: `String(QDir(...).absolutePath())`)

Note: `MaRaClusterAdapter.cpp` uses `.toQString()` for String interop but its Qt include (`<QDir>`) is stale — handle in both Section 5 and Section 6.

For tools that don't link openms_gui, inline the conversions:
- `str.toQString()` → `QString::fromStdString(str)`
- `String(some_qstring)` → `some_qstring.toStdString()` (assigns to `std::string`/`String`)

### 6. CMake: Explicit Qt Linking for TOPP Tools and Tests

Removing `Qt6::Core` from libOpenMS's PUBLIC dependencies means TOPP tools and tests that directly `#include` Qt headers will lose their transitive Qt dependency and fail to compile.

**Fix:** Add `Qt6::Core` as an explicit link dependency to each affected target individually. For targets with only stale includes, remove the include instead.

**Stale Qt includes (remove the include, no Qt link needed):**
- TOPP: `PSMFeatureExtractor.cpp` (QFile, QDir, QProcess), `NucleicAcidSearchEngine.cpp` (QProcess), `MapStatistics.cpp` (QString), `FeatureFinderMultiplex.cpp` (QDir), `LuciphorAdapter.cpp` (QProcessEnvironment), `MaRaClusterAdapter.cpp` (QDir), `QCExporter.cpp`, `QCImporter.cpp`, `QCMerger.cpp`, `Resampler.cpp`
- Tests: `MzMLSqliteSwathHandler_test` (QFile), `StringListUtils_test` (QStringList), `SqMassFile_test` (QFile), `OPXLHelper_test` (QStringList), `DataValue_test` (QString), `SiriusFragmentAnnotation_test` (QDir, QString), `XQuestResultXMLFile_test` (QStringList)

**Affected TOPP tools (link only OpenMS, actively use Qt — need explicit `Qt6::Core`):**
- `QCEmbedder`, `QCExtractor`, `QCShrinker` — QByteArray, QFile, QFileInfo, QString
- `OpenSwathWorkflow` — QDir
- Plus tools from Section 5 that use Qt beyond String interop: `MetaProSIP`, `GenericWrapper`, `IDRipper`, `MzMLSplitter`, `CometAdapter`, `OpenSwathFileSplitter`, `AssayGeneratorMetaboSirius`

**Affected tests (link only OpenMS, actively use Qt — need explicit `Qt6::Core`):**
- `PeptideIndexing_test` (QStringList), `FIAMSScheduler_test` (QDir)
- Plus the 3 tests from Section 7 that use String/Qt interop: `MzMLSqliteHandler_test`, `PythonInfo_test`, `ToolDescriptionFile_test`

### 7. Test Changes (`src/tests/`)

- **`String_test.cpp`** — Remove test cases for `String(const QString&)` and `toQString()`. All other tests untouched.
- **`MzMLSqliteHandler_test.cpp`**, **`PythonInfo_test.cpp`**, **`ToolDescriptionFile_test.cpp`** — Inline conversions directly (`QString::fromStdString()` / `.toStdString()`) since these tests don't link openms_gui.

## Result

- `libOpenMS` has zero Qt headers, zero Qt symbols, zero Qt link dependency
- CMake configuration no longer requires Qt for core-lib-only consumers
- pyOpenMS builds without Qt
- GUI code continues working via `toQString()`/`fromQString()` free functions in `Qt5Port.h`
- TOPP tools and tests use inlined conversions (`QString::fromStdString()` / `.toStdString()`)
- Breaking change for downstream C++ code that uses `String(qstr)` or `str.toQString()` — must switch to direct Qt string conversions or `Qt5Port.h` free functions. Document in CHANGELOG.

## File Impact Summary

| Area | Files | Change Type |
|------|-------|-------------|
| Core lib | 4 | Remove Qt API + dependency |
| CMake config | 3 | Make Qt optional for core lib |
| Compat layer (Qt5Port.h) | 1 | Add free functions + explicit `<QString>` include |
| GUI lib | ~60 | Mechanical migration (free functions, includes headers) |
| TOPP tools | 12 | Mechanical migration (inline conversions or free functions) |
| Stale Qt includes | ~17 | Remove unused Qt includes (TOPP + tests) |
| CMake (TOPP) | ~11 | Add explicit `Qt6::Core` link dependency |
| CMake (tests) | ~5 | Add explicit `Qt6::Core` link dependency |
| Tests | 4 | Remove Qt tests + inline conversions |
