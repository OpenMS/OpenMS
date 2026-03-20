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

### 2. Compatibility Layer (`src/openms_gui/`)

**`src/openms_gui/include/OpenMS/VISUAL/MISC/Qt5Port.h`** — Add two free functions:

```cpp
/// Convert OpenMS::String to QString
inline QString toQString(const String& s) { return QString::fromStdString(s); }

/// Construct OpenMS::String from QString
inline String fromQString(const QString& s) { return s.toStdString(); }
```

This header already includes `<OpenMS/DATASTRUCTURES/String.h>` and `<QStringList>`.

### 3. GUI Library Migration (`src/openms_gui/`)

~60 source files that call `.toQString()` or use `String(qstr)` implicit conversion:
- Replace `str.toQString()` with `toQString(str)` (free function from `Qt5Port.h`)
- Replace `String(some_qstring)` with `fromQString(some_qstring)` (free function from `Qt5Port.h`)
- Add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>` where not already included

These files already link `OpenMS_GUI` and have the include path available.

### 4. TOPP Tool Migration (`src/topp/`)

11 files that call `.toQString()` or use `String(qstr)`:

**Link openms_gui (can use `Qt5Port.h`):** `ExecutePipeline.cpp`, `ImageCreator.cpp`, `INIUpdater.cpp`

**Link only OpenMS (inline conversions directly):** `MetaProSIP.cpp`, `GenericWrapper.cpp`, `IDRipper.cpp`, `MzMLSplitter.cpp`, `CometAdapter.cpp`, `MaRaClusterAdapter.cpp`, `OpenSwathFileSplitter.cpp`, `AssayGeneratorMetaboSirius.cpp`

For tools that don't link openms_gui, inline the conversions:
- `str.toQString()` → `QString::fromStdString(str)`
- `String(some_qstring)` → `some_qstring.toStdString()` (assigns to `std::string`/`String`)

### 5. CMake: Explicit Qt Linking for TOPP Tools and Tests

Removing `Qt6::Core` from libOpenMS's PUBLIC dependencies means TOPP tools and tests that directly `#include` Qt headers (QDir, QFile, QByteArray, QString, etc.) will lose their transitive Qt dependency and fail to compile.

**Fix:** Add `Qt6::Core` as an explicit link dependency to each affected target individually in `src/topp/CMakeLists.txt` and the relevant test CMakeLists.

**Stale Qt includes (remove the include, no Qt link needed):**
- TOPP: `PSMFeatureExtractor.cpp` (QFile, QDir, QProcess — unused), `NucleicAcidSearchEngine.cpp` (QProcess — unused)
- Tests: `MzMLSqliteSwathHandler_test` (QFile), `StringListUtils_test` (QStringList), `SqMassFile_test` (QFile), `OPXLHelper_test` (QStringList), `DataValue_test` (QString), `SiriusFragmentAnnotation_test` (QDir, QString), `XQuestResultXMLFile_test` (QStringList) — all unused

**Affected TOPP tools (link only OpenMS, actively use Qt — need explicit `Qt6::Core`):**
- `QCEmbedder`, `QCExporter`, `QCExtractor`, `QCImporter`, `QCMerger`, `QCShrinker` — QByteArray, QFile, QFileInfo, QString
- `MapStatistics` — QtCore/QString
- `FeatureFinderMultiplex` — QDir
- `LuciphorAdapter` — QProcessEnvironment
- `OpenSwathWorkflow` — QDir
- Plus the 8 tools from Section 4 that use Qt beyond String interop: `MetaProSIP`, `GenericWrapper`, `IDRipper`, `MzMLSplitter`, `CometAdapter`, `MaRaClusterAdapter`, `OpenSwathFileSplitter`, `AssayGeneratorMetaboSirius`

**Affected tests (link only OpenMS, actively use Qt — need explicit `Qt6::Core`):**
- `PeptideIndexing_test` (QStringList), `FIAMSScheduler_test` (QDir)
- Plus the 3 tests from Section 6 that use String/Qt interop: `MzMLSqliteHandler_test`, `PythonInfo_test`, `ToolDescriptionFile_test`

### 6. Test Changes (`src/tests/`)

- **`String_test.cpp`** — Remove test cases for `String(const QString&)` and `toQString()`. All other tests untouched.
- **`MzMLSqliteHandler_test.cpp`**, **`PythonInfo_test.cpp`**, **`ToolDescriptionFile_test.cpp`** — Inline conversions directly (`QString::fromStdString()` / `.toStdString()`) since these tests don't link openms_gui.

## Result

- `libOpenMS` has zero Qt headers, zero Qt symbols, zero Qt link dependency
- pyOpenMS builds without Qt
- GUI code continues working via `toQString()`/`fromQString()` free functions in `Qt5Port.h`
- TOPP tools and tests use inlined conversions (`QString::fromStdString()` / `.toStdString()`)
- Breaking change for downstream C++ code that uses `String(qstr)` or `str.toQString()` — must switch to direct Qt string conversions or `Qt5Port.h` free functions. Document in CHANGELOG.

## File Impact Summary

| Area | Files | Change Type |
|------|-------|-------------|
| Core lib | 4 | Remove Qt API + dependency |
| Compat layer (Qt5Port.h) | 1 | Add free functions |
| GUI lib | ~60 | Mechanical migration (free functions) |
| TOPP tools | 11 | Mechanical migration (inline conversions or free functions) |
| Stale Qt includes | ~11 | Remove unused Qt includes |
| CMake (TOPP) | ~18 | Add explicit `Qt6::Core` link dependency |
| CMake (tests) | ~5 | Add explicit `Qt6::Core` link dependency |
| Tests | 4 | Remove Qt tests + inline conversions |
