# Design: Remove Qt from TOPP Tools

## Goal

Make `src/topp/` 100% Qt-free. All Qt usage moves to `src/openms_gui/` or is replaced with `std::filesystem`, `std::regex`, and existing OpenMS core utilities (`File`, `ExternalProcess`, `FileHandler::stripExtension`, `Base64`).

## Current State

- OpenMS core library (`src/openms/`) is already Qt-free (commit `c23b814893`).
- 14 of 148 TOPP tools still use Qt:
  - **10 tools** link `Qt6::Core` for file/path/regex/process operations
  - **4 tools** link `OpenMS_GUI` (ExecutePipeline, INIUpdater, ImageCreator, Resampler)
- Resampler's GUI dependency is a dead include — it doesn't actually use any GUI code.

## Phase 1: Replace Qt with std/OpenMS utilities (11 tools)

### 1.1 Resampler — remove dead include

- Remove `#include <OpenMS/VISUAL/MultiGradient.h>` (unused)
- Move Resampler from `TOPP_executables_with_GUIlib` to `TOPP_executables` in `src/topp/executables.cmake`
- Remove Resampler from the GUI-only erasure list in `src/openms/source/APPLICATIONS/ToolHandler.cpp` (lines ~196-207)
- Remove the `if(WITH_GUI)` gate around Resampler tests in `src/tests/topp/CMakeLists.txt` (lines ~2021-2029)

### 1.2 Path operations — `std::filesystem` + `OpenMS::File` (7 tools)

**Replacement mapping:**

| Qt API | Replacement |
|--------|-------------|
| `QDir(path).absolutePath()` | `File::absolutePath(path)` |
| `QFileInfo(f).baseName()` | `FileHandler::stripExtension(File::basename(f))` — handles compound extensions like `.mzML.gz` correctly |
| `QFileInfo(f).completeBaseName()` | `std::filesystem::path(f).stem().string()` — strips only last extension, matching `completeBaseName()` semantics |
| `QFileInfo(f).absoluteFilePath()` | `File::absolutePath(f)` |
| `QDir::toNativeSeparators(p)` | `std::filesystem::path(p).make_preferred().string()` |
| `QFile(f).size()` | `std::filesystem::file_size(f)` |
| `QFile::readAll()` + `QByteArray::toBase64()` | `std::ifstream` binary read + `Base64::encode()` from `OpenMS/FORMAT/Base64.h` |

**Important:** Do NOT use `std::filesystem::path::stem()` as a replacement for `QFileInfo::baseName()`. They differ on multi-dot filenames (e.g., `sample.mzML.gz` — `baseName()` returns `sample`, `stem()` returns `sample.mzML`). Use `FileHandler::stripExtension(File::basename(path))` instead, which understands OpenMS file types.

**Tool-specific changes:**

- **MzMLSplitter**: `QFile(f).size()` → `std::filesystem::file_size(f)`. Remove `<QFile>` include.
- **IDRipper**: `QFileInfo::absoluteFilePath()` → `File::absolutePath()`. `QFileInfo::completeBaseName()` → `std::filesystem::path::stem()` (safe here — `completeBaseName` and `stem` have same semantics). `QDir::toNativeSeparators()` → `std::filesystem::path::make_preferred()`. Remove `<QDir>` include.
- **OpenSwathFileSplitter**: `QDir::absolutePath()` → `File::absolutePath()`. `QFileInfo::baseName()` → `FileHandler::stripExtension(File::basename(f))`. Remove `<QDir>` include. Do NOT use native path separators — code relies on forward slashes.
- **OpenSwathWorkflow**: `QDir::absolutePath()` → `File::absolutePath()`. Remove `<QDir>` include. Do NOT use native path separators.
- **QCEmbedder**: `QFileInfo::baseName()` → `FileHandler::stripExtension(File::basename(f))`. `QFile::readAll()` → `std::ifstream` binary read. `QByteArray::toBase64()` → `Base64::encode()`. Remove Qt includes.
- **QCExtractor**: `QFileInfo::baseName()` → `FileHandler::stripExtension(File::basename(f))`. Remove Qt includes.
- **QCShrinker**: `QFileInfo::baseName()` → `FileHandler::stripExtension(File::basename(f))`. Remove Qt includes.

### 1.3 Directory iteration — AssayGeneratorMetaboSirius

- Replace `QDirIterator` with `std::filesystem::directory_iterator`.
- Replace `QString::fromStdString()` / `.toStdString()` with direct `std::string` / `std::filesystem::path` usage.
- Add `std::error_code` overload or try/catch for `directory_iterator` construction, since it throws on invalid paths where `QDirIterator` silently returns empty. Match existing behavior (treat bad path as empty).

### 1.4 Regex — CometAdapter

- Replace `QRegularExpression` with `std::regex`.
- Replace `QString::fromStdString().remove(QRegularExpression(...))` with `std::regex_replace()`.
- The patterns used are simple enough that no syntax differences apply.

### 1.5 Process execution — MetaProSIP

This is the broadest single-tool refactor. Full Qt surface used:

| Qt API | Replacement |
|--------|-------------|
| `QProcess` (start, waitForFinished, exitCode, error, exitStatus) | `ExternalProcess::run()` — returns `RETURNSTATE` enum covering `SUCCESS`, `NONZERO_EXIT`, `CRASH`, `FAILED_TO_START` |
| `QProcess::systemEnvironment()` + `setEnvironment()` | `ExternalProcess` `env` parameter (`std::map<String, String>`) |
| `QProcess::MergedChannels` + `readAllStandardOutput()` | Callback-based: pass lambda that accumulates into `std::string` |
| `QFile::copy()` | `File::copy()` |
| `QFile::remove()` | `File::remove()` |
| `QDir::absolutePath()` | `File::absolutePath()` |
| `QDir::exists()` / `QDir::mkpath()` | `File::exists()` / `std::filesystem::create_directories()` |
| `QFileInfo::baseName()` | `FileHandler::stripExtension(File::basename(f))` |
| `QString` typed parameters | `std::string` / `OpenMS::String` |

### 1.6 CMake changes (Phase 1)

- **`src/topp/executables.cmake`**: Move Resampler from `TOPP_executables_with_GUIlib` to `TOPP_executables`.
- **`src/topp/CMakeLists.txt`**: Remove the entire `if(TARGET Qt6::Core)` block that links 10 tools to `Qt6::Core`.
- **`src/openms/source/APPLICATIONS/ToolHandler.cpp`**: Remove Resampler from the `GUI_tools` erasure list.
- **`src/tests/topp/CMakeLists.txt`**: Remove `if(WITH_GUI)` gate around Resampler tests.

## Phase 2: Move 3 GUI tools to `src/openms_gui/`

### 2.1 File moves

- `src/topp/ExecutePipeline.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/ExecutePipeline.cpp`
- `src/topp/INIUpdater.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/INIUpdater.cpp`
- `src/topp/ImageCreator.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/ImageCreator.cpp`

### 2.2 CMake changes

- **`src/topp/executables.cmake`**: Remove ExecutePipeline, INIUpdater, ImageCreator from `TOPP_executables_with_GUIlib`. Remove the `TOPP_executables_with_GUIlib` variable entirely (it will be empty).
- **`src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/executables.cmake`**: Add ExecutePipeline, INIUpdater, ImageCreator to `GUI_executables` list.
- **`src/openms_gui/CMakeLists.txt`**: These 3 tools need to link both `OpenMS` and `OpenMS_GUI`. The existing GUI executable build loop already links the Visual library; verify that `OpenMS` is also linked (it should be transitively via `OpenMS_GUI`).
- **`src/topp/CMakeLists.txt`**: Remove the entire `if(WITH_GUI)` block for `TOPP_executables_with_GUIlib` (no longer needed).

### 2.3 Tool registration

- **`src/openms/source/APPLICATIONS/ToolHandler.cpp`**: The 3 moved tools must remain in the TOPP tool registry so they appear in `-list_utils` output and are recognized by `TOPPBase`. Currently they're in `TOPP_executables_with_GUIlib` which gets merged into the tool map. After the move, they need to be registered as GUI-only tools (present in the map but erased when `WITH_GUI` is off — the existing `#ifndef WITH_GUI` erasure pattern).

### 2.4 Test and automation impact

- **`TOPP_TOOLS` variable**: Currently exported from `src/topp/CMakeLists.txt` and drives:
  - Auto-generated `-write_ini` / `-write_ctd` tests in `src/tests/topp/CMakeLists.txt`
  - CWL generation in `cmake/cwl_generation.cmake`
- The 3 moved tools must be added to this variable (or a parallel `GUI_TOPP_TOOLS` variable) so their INI/CTD/CWL artifacts are still generated when `WITH_GUI=ON`.
- Existing tool-specific tests in `src/tests/topp/CMakeLists.txt` reference tools by binary name via `${TOPP_BIN_PATH}/ToolName`, not source path. As long as the built binaries end up in the same bin directory, tests work unchanged. Keep existing `if(WITH_GUI)` gates on these tests.

### 2.5 Qt find_package cleanup

- **`cmake/cmake_findExternalLibs.cmake`**: The `find_package(Qt6 Core QUIET)` on line ~342 is currently outside the `if(WITH_GUI)` block because TOPP tools needed it. After Phase 1+2, move it inside the `if(WITH_GUI)` block so Qt is only searched when building GUI.
- **`cmake/knime_package_support.cmake`**: Uses `Qt6::qmake` to query Qt install paths (~line 68). This is only relevant for KNIME packaging which implies GUI. Guard with `if(WITH_GUI)` or `if(TARGET Qt6::qmake)`.

## Verification

- Build with `WITH_GUI=ON`: all tools compile, `ctest` passes
- Build with `WITH_GUI=OFF`: all non-GUI TOPP tools compile without Qt, GUI tools are skipped
- `grep -r "Qt" src/topp/` returns no hits
- Verify CWL generation still produces artifacts for all tools (including the 3 moved ones)
- MetaProSIP: test R script invocation error paths (failed to start, crash, non-zero exit)

## Risk Assessment

- **MetaProSIP**: Highest risk — broadest refactor, callback-based output capture is a design change. Needs careful testing.
- **baseName semantics**: Medium risk — using `FileHandler::stripExtension(File::basename())` is safer than `stem()` but must be verified against actual file inputs in tests.
- **CometAdapter regex**: Low risk — simple pattern, `std::regex` handles it identically.
- **Tool registration after move**: Medium risk — must ensure ToolHandler, CWL, and INI/CTD generation all still find the 3 moved tools.
