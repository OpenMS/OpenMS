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

Minimal surgical removal: remove Qt from the core lib's public API and provide free-function replacements in `openms_gui`'s existing `Qt5Port.h` compatibility header.

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
- `#include <QtCore/QDateTime>` (line 46, stale/unused)

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

~40 source files that call `.toQString()` or use `String(qstr)` implicit conversion:
- Replace `str.toQString()` with `toQString(str)`
- Replace `String(some_qstring)` with `fromQString(some_qstring)`
- Add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>` where not already included

### 4. TOPP Tool Migration (`src/topp/`)

11 files that call `.toQString()` or use `String(qstr)`:
- `MetaProSIP.cpp` (~30 calls)
- `GenericWrapper.cpp` (~15 calls)
- `INIUpdater.cpp` (~10 calls)
- `ExecutePipeline.cpp` (~7 calls)
- `IDRipper.cpp`
- `ImageCreator.cpp`
- `MzMLSplitter.cpp`
- `CometAdapter.cpp`
- `MaRaClusterAdapter.cpp`
- `OpenSwathFileSplitter.cpp`
- `AssayGeneratorMetaboSirius.cpp`

Same migration pattern: add `#include <OpenMS/VISUAL/MISC/Qt5Port.h>`, replace member calls with free functions.

### 5. Test Changes (`src/tests/`)

- **`String_test.cpp`** — Remove test cases for `String(const QString&)` and `toQString()`. All other tests untouched.
- **`MzMLSqliteHandler_test.cpp`**, **`PythonInfo_test.cpp`**, **`ToolDescriptionFile_test.cpp`** — Replace `.toQString()` calls with `toQString()` free function via `Qt5Port.h`.

## Result

- `libOpenMS` has zero Qt headers, zero Qt symbols, zero Qt link dependency
- pyOpenMS builds without Qt
- GUI code and TOPP tools continue working via `toQString()`/`fromQString()` free functions in `Qt5Port.h`
- Breaking change for downstream C++ code that uses `String(qstr)` or `str.toQString()` — must switch to `Qt5Port.h` free functions

## File Impact Summary

| Area | Files | Change Type |
|------|-------|-------------|
| Core lib | 4 | Remove Qt API + dependency |
| Compat layer (Qt5Port.h) | 1 | Add free functions |
| GUI lib | ~40 | Mechanical migration |
| TOPP tools | 11 | Mechanical migration |
| Tests | 4 | Remove Qt tests + migrate callers |
