# Remove Qt from TOPP Tools — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `src/topp/` 100% Qt-free by replacing Qt usage with `std::filesystem`, `std::regex`, and existing OpenMS utilities, then moving 3 GUI-entangled tools into `src/openms_gui/`.

**Architecture:** Two-phase approach. Phase 1 replaces Qt6::Core usage in 11 TOPP tools with standard C++ and existing OpenMS APIs (`File`, `ExternalProcess`, `FileHandler::stripExtension`). Phase 2 moves 3 genuinely GUI-dependent tools (ExecutePipeline, INIUpdater, ImageCreator) into the GUI library directory and cleans up CMake so Qt is only found when `WITH_GUI=ON`.

**Tech Stack:** C++20, std::filesystem, std::regex, OpenMS::File, OpenMS::ExternalProcess, OpenMS::FileHandler, OpenMS::Base64, CMake

**Spec:** `docs/superpowers/specs/2026-03-30-remove-qt-from-topp-tools-design.md`

---

## File Map

### Phase 1 — Files modified:
- `src/topp/Resampler.cpp` — remove dead include
- `src/topp/MzMLSplitter.cpp` — replace QFile::size()
- `src/topp/IDRipper.cpp` — replace QFileInfo/QDir path ops
- `src/topp/OpenSwathFileSplitter.cpp` — replace QDir/QFileInfo
- `src/topp/OpenSwathWorkflow.cpp` — replace QDir
- `src/topp/QCExtractor.cpp` — replace QFileInfo::baseName()
- `src/topp/QCShrinker.cpp` — replace QFileInfo::baseName()
- `src/topp/QCEmbedder.cpp` — replace QFileInfo, QFile, QByteArray
- `src/topp/AssayGeneratorMetaboSirius.cpp` — replace QDirIterator
- `src/topp/CometAdapter.cpp` — replace QRegularExpression
- `src/topp/MetaProSIP.cpp` — replace QProcess, QFile, QDir, QFileInfo, QString
- `src/topp/executables.cmake` — move Resampler to regular list, remove GUI list
- `src/topp/CMakeLists.txt` — remove Qt6::Core linking block, remove WITH_GUI block
- `src/openms/source/APPLICATIONS/ToolHandler.cpp` — update GUI_tools erasure list
- `src/tests/topp/CMakeLists.txt` — ungate Resampler tests

### Phase 2 — Files moved/modified:
- `src/topp/ExecutePipeline.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/ExecutePipeline.cpp`
- `src/topp/INIUpdater.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/INIUpdater.cpp`
- `src/topp/ImageCreator.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/ImageCreator.cpp`
- `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/executables.cmake` — add 3 tools
- `src/openms_gui/CMakeLists.txt` — verify linking
- `cmake/cmake_findExternalLibs.cmake` — move Qt6::Core find inside WITH_GUI
- `cmake/knime_package_support.cmake` — guard Qt6::qmake usage

---

### Task 1: Resampler — remove dead include and ungate from GUI

**Files:**
- Modify: `src/topp/Resampler.cpp:15`
- Modify: `src/topp/executables.cmake:155-168`
- Modify: `src/openms/source/APPLICATIONS/ToolHandler.cpp:196-207`
- Modify: `src/tests/topp/CMakeLists.txt:2021-2029`

- [ ] **Step 1: Remove dead include from Resampler.cpp**

In `src/topp/Resampler.cpp`, remove line 15:

```cpp
// REMOVE this line:
#include <OpenMS/VISUAL/MultiGradient.h>
```

- [ ] **Step 2: Move Resampler from GUI list to regular list in executables.cmake**

In `src/topp/executables.cmake`, add `Resampler` to the `TOPP_executables` list (alphabetically, around line 120 near other R-tools), and remove it from `TOPP_executables_with_GUIlib` (line 164):

```cmake
# In TOPP_executables list, add after ProteomicsLFQ (or wherever alphabetical):
Resampler
```

```cmake
# TOPP_executables_with_GUIlib becomes:
set(TOPP_executables_with_GUIlib
ExecutePipeline
# util category
ImageCreator
INIUpdater
)
```

- [ ] **Step 3: Remove Resampler from GUI_tools erasure in ToolHandler.cpp**

In `src/openms/source/APPLICATIONS/ToolHandler.cpp`, remove `"Resampler",` from the `GUI_tools` list (line 202):

```cpp
#ifndef WITH_GUI
    StringList GUI_tools = {
      "ExecutePipeline",
      "ImageCreator",
      "INIUpdater",
    };
```

- [ ] **Step 4: Ungate Resampler tests**

In `src/tests/topp/CMakeLists.txt`, remove the `if(WITH_GUI)` / `endif()` wrapper around Resampler tests (lines 2021 and 2029), keeping the test definitions:

```cmake
#------------------------------------------------------------------------------
# Resampler tests
  add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out ${TESTS_TEMP_DIR}/Resampler.mzML -sampling_rate 0.3)
  add_test("TOPP_Resampler_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 ${TESTS_TEMP_DIR}/Resampler.mzML -in2 ${DATA_DIR_TOPP}/Resampler_1_output.mzML )
  set_tests_properties("TOPP_Resampler_1_out1" PROPERTIES DEPENDS "TOPP_Resampler_1")

  add_test("TOPP_Resampler_2" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out ${TESTS_TEMP_DIR}/Resampler2_ppm.mzML -ppm -sampling_rate 2000 -min_int_cutoff 3900)
  add_test("TOPP_Resampler_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 ${TESTS_TEMP_DIR}/Resampler2_ppm.mzML -in2 ${DATA_DIR_TOPP}/Resampler_2_output.mzML )
  set_tests_properties("TOPP_Resampler_2_out1" PROPERTIES DEPENDS "TOPP_Resampler_2")
```

- [ ] **Step 5: Build and verify**

Run: `cmake --build OpenMS-build --target Resampler -j$(nproc)`
Expected: Compiles without Qt, no MultiGradient.h dependency.

- [ ] **Step 6: Commit**

```bash
git add src/topp/Resampler.cpp src/topp/executables.cmake src/openms/source/APPLICATIONS/ToolHandler.cpp src/tests/topp/CMakeLists.txt
git commit -m "refactor(Resampler): remove dead GUI include, ungate from WITH_GUI"
```

---

### Task 2: MzMLSplitter — replace QFile::size()

**Files:**
- Modify: `src/topp/MzMLSplitter.cpp:13,98-100`

- [ ] **Step 1: Replace Qt include and usage**

In `src/topp/MzMLSplitter.cpp`:

Remove line 13:
```cpp
// REMOVE:
#include <QFile>
```

Add filesystem include (near other includes):
```cpp
#include <filesystem>
```

Replace lines 98-100:
```cpp
// OLD:
      QFile mzml_file(QString::fromStdString(in));
      // use float here to avoid too many decimals in output below:
      float total_size = mzml_file.size();

// NEW:
      // use float here to avoid too many decimals in output below:
      float total_size = static_cast<float>(std::filesystem::file_size(in));
```

- [ ] **Step 2: Build and verify**

Run: `cmake --build OpenMS-build --target MzMLSplitter -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 3: Commit**

```bash
git add src/topp/MzMLSplitter.cpp
git commit -m "refactor(MzMLSplitter): replace QFile::size() with std::filesystem::file_size()"
```

---

### Task 3: QCExtractor and QCShrinker — replace QFileInfo::baseName()

**Files:**
- Modify: `src/topp/QCExtractor.cpp:15-18,112`
- Modify: `src/topp/QCShrinker.cpp:15-18,110`

- [ ] **Step 1: Replace Qt in QCExtractor.cpp**

Remove lines 15-18:
```cpp
// REMOVE all four:
#include <QByteArray>
#include <QFile>
#include <QString>
#include <QFileInfo>
```

Add includes:
```cpp
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
```

Replace line 112:
```cpp
// OLD:
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName().toStdString();

// NEW:
      target_run = FileHandler::stripExtension(File::basename(target_file));
```

- [ ] **Step 2: Replace Qt in QCShrinker.cpp**

Same pattern. Remove lines 15-18 (same Qt includes). Add same OpenMS includes.

Replace line 110:
```cpp
// OLD:
      target_run = QFileInfo(QString::fromStdString(target_file)).baseName().toStdString();

// NEW:
      target_run = FileHandler::stripExtension(File::basename(target_file));
```

- [ ] **Step 3: Build and verify**

Run: `cmake --build OpenMS-build --target QCExtractor QCShrinker -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 4: Commit**

```bash
git add src/topp/QCExtractor.cpp src/topp/QCShrinker.cpp
git commit -m "refactor(QCExtractor,QCShrinker): replace QFileInfo::baseName() with FileHandler::stripExtension()"
```

---

### Task 4: QCEmbedder — replace QFileInfo, QFile, QByteArray

**Files:**
- Modify: `src/topp/QCEmbedder.cpp:19-22,137,162-168`

- [ ] **Step 1: Replace includes**

Remove lines 19-22:
```cpp
// REMOVE:
#include <QByteArray>
#include <QFile>
#include <QString>
#include <QFileInfo>
```

Add:
```cpp
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <sstream>
#include <OpenMS/FORMAT/Base64.h>
```

- [ ] **Step 2: Replace baseName usage (line 137)**

```cpp
// OLD:
    target_run = QFileInfo(QString::fromStdString(target_file)).baseName().toStdString();

// NEW:
    target_run = FileHandler::stripExtension(File::basename(target_file));
```

- [ ] **Step 3: Replace file read + base64 encoding (lines 162-168)**

```cpp
// OLD:
  QFile f(plot_file.c_str());
  String plot_b64;
  if (f.open(QIODevice::ReadOnly))
  {
    QByteArray ba = f.readAll();
    f.close();
    plot_b64 = String(ba.toBase64().toStdString());
  }

// NEW:
  String plot_b64;
  {
    std::ifstream f(plot_file, std::ios::binary);
    if (f.is_open())
    {
      std::ostringstream oss;
      oss << f.rdbuf();
      std::string raw = oss.str();
      std::vector<unsigned char> data(raw.begin(), raw.end());
      Base64::encode(data, Base64::BYTEORDER_LITTLEENDIAN, plot_b64, false);
    }
  }
```

Note: Check that `Base64::encode` with `unsigned char` and `BYTEORDER_LITTLEENDIAN` produces standard base64 matching `QByteArray::toBase64()`. If the OpenMS `Base64` class uses a different encoding scheme (it may be tailored for numerical arrays), a simple standalone base64 encoder may be needed instead. Verify by comparing output on a small test PNG.

- [ ] **Step 4: Build and verify**

Run: `cmake --build OpenMS-build --target QCEmbedder -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 5: Commit**

```bash
git add src/topp/QCEmbedder.cpp
git commit -m "refactor(QCEmbedder): replace QFile/QByteArray with std::ifstream and Base64"
```

---

### Task 5: IDRipper — replace QFileInfo/QDir path operations

**Files:**
- Modify: `src/topp/IDRipper.cpp:16,102,133-134,140,144,147,150`

- [ ] **Step 1: Replace include**

Remove line 16:
```cpp
// REMOVE:
#include <QDir>
```

Add:
```cpp
#include <OpenMS/SYSTEM/File.h>
#include <filesystem>
```

- [ ] **Step 2: Replace absoluteFilePath (line 102)**

```cpp
// OLD:
    String output_directory = QFileInfo(QString::fromStdString(out_dir)).absoluteFilePath().toStdString();

// NEW:
    String output_directory = File::absolutePath(out_dir);
```

- [ ] **Step 3: Replace completeBaseName and path construction (lines 133-150)**

```cpp
// OLD:
      QString output = QString::fromStdString(output_directory);

      String out_fname;
      if (numeric_filenames)
      {
        String s_ident_run_idx = split_ident_runs ? '_' + String(rfi.ident_run_idx) : "";
        String s_file_origin_idx = '_' + String(rfi.file_origin_idx);
        out_fname = QFileInfo(QString::fromStdString(file_name)).completeBaseName().toStdString() + s_ident_run_idx + s_file_origin_idx + ".idXML";
      }
      else
      {
        out_fname = QFileInfo(QString::fromStdString(rfi.out_basename)).completeBaseName().toStdString() + ".idXML";
      }

      String out = QDir::toNativeSeparators(output.append(QString("/")).append(QString::fromStdString(out_fname))).toStdString();
      OPENMS_LOG_INFO << "Storing file: '" << out << "'." << std::endl;

      QDir dir(QString::fromStdString(output_directory));

// NEW:
      String out_fname;
      if (numeric_filenames)
      {
        String s_ident_run_idx = split_ident_runs ? '_' + String(rfi.ident_run_idx) : "";
        String s_file_origin_idx = '_' + String(rfi.file_origin_idx);
        out_fname = std::filesystem::path(std::string(file_name)).stem().string() + s_ident_run_idx + s_file_origin_idx + ".idXML";
      }
      else
      {
        out_fname = std::filesystem::path(std::string(rfi.out_basename)).stem().string() + ".idXML";
      }

      String out = (std::filesystem::path(std::string(output_directory)) / std::string(out_fname)).make_preferred().string();
      OPENMS_LOG_INFO << "Storing file: '" << out << "'." << std::endl;
```

Note: `std::filesystem::path::stem()` is correct here because IDRipper uses `completeBaseName()` (strips only last extension), which matches `stem()` semantics.

- [ ] **Step 4: Build and verify**

Run: `cmake --build OpenMS-build --target IDRipper -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 5: Commit**

```bash
git add src/topp/IDRipper.cpp
git commit -m "refactor(IDRipper): replace QDir/QFileInfo with std::filesystem and File::absolutePath()"
```

---

### Task 6: OpenSwathFileSplitter and OpenSwathWorkflow — replace QDir

**Files:**
- Modify: `src/topp/OpenSwathFileSplitter.cpp:18,97,99-100`
- Modify: `src/topp/OpenSwathWorkflow.cpp:71,638`

- [ ] **Step 1: Replace Qt in OpenSwathFileSplitter.cpp**

Remove line 18:
```cpp
// REMOVE:
#include <QDir>
```

Add:
```cpp
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
```

Replace lines 97-100:
```cpp
// OLD:
    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // (do not use QDir::separator(), since its platform specific (/ or \) while absolutePath() will always use '/')
    String tmp_dir = String(QDir(getStringOption_("outputDirectory").c_str()).absolutePath().toStdString()).ensureLastChar('/');

    QFileInfo fi(QString::fromStdString(file_in));
    String tmp = tmp_dir + String(fi.baseName().toStdString());

// NEW:
    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // File::absolutePath() always uses '/' separators
    String tmp_dir = File::absolutePath(getStringOption_("outputDirectory")).ensureLastChar('/');

    String tmp = tmp_dir + FileHandler::stripExtension(File::basename(file_in));
```

- [ ] **Step 2: Replace Qt in OpenSwathWorkflow.cpp**

Remove line 71:
```cpp
// REMOVE:
#include <QDir>
```

Add:
```cpp
#include <OpenMS/SYSTEM/File.h>
```

Replace line 638:
```cpp
// OLD:
    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // (do not use QDir::separator(), since its platform specific (/ or \) while absolutePath() will always use '/')
    String tmp_dir = String(QDir(getStringOption_("tempDirectory").c_str()).absolutePath().toStdString()).ensureLastChar('/');

// NEW:
    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // File::absolutePath() always uses '/' separators
    String tmp_dir = File::absolutePath(getStringOption_("tempDirectory")).ensureLastChar('/');
```

- [ ] **Step 3: Build and verify**

Run: `cmake --build OpenMS-build --target OpenSwathFileSplitter OpenSwathWorkflow -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 4: Commit**

```bash
git add src/topp/OpenSwathFileSplitter.cpp src/topp/OpenSwathWorkflow.cpp
git commit -m "refactor(OpenSwath*): replace QDir with File::absolutePath()"
```

---

### Task 7: AssayGeneratorMetaboSirius — replace QDirIterator

**Files:**
- Modify: `src/topp/AssayGeneratorMetaboSirius.cpp:26-28,195-198`

- [ ] **Step 1: Replace includes**

Remove lines 26-28:
```cpp
// REMOVE:
#include <QtCore/QDir>
#include <QtCore/QDirIterator>
#include <QtCore/QString>
```

Add:
```cpp
#include <filesystem>
```

- [ ] **Step 2: Replace QDirIterator usage (lines 192-198)**

```cpp
// OLD:
    // Get all subdirectories within the SIRIUS project directory
    //-------------------------------------------------------------
    std::vector<String> subdirs;
    QDirIterator it(QString::fromStdString(sirius_project_directory), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
    while (it.hasNext())
    {
      subdirs.emplace_back(it.next().toStdString());

// NEW:
    // Get all subdirectories within the SIRIUS project directory
    //-------------------------------------------------------------
    std::vector<String> subdirs;
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(std::string(sirius_project_directory), ec))
    {
      if (entry.is_directory())
      {
        subdirs.emplace_back(entry.path().string());
      }
```

Note: Using the `std::error_code` overload so that invalid paths produce an empty iteration (matching QDirIterator behavior) instead of throwing.

Also update the closing brace — the old `while` loop had a single closing `}`. The new `for` loop also has a single closing `}` after the `if` block, so check surrounding code for correct brace matching.

- [ ] **Step 3: Build and verify**

Run: `cmake --build OpenMS-build --target AssayGeneratorMetaboSirius -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 4: Commit**

```bash
git add src/topp/AssayGeneratorMetaboSirius.cpp
git commit -m "refactor(AssayGeneratorMetaboSirius): replace QDirIterator with std::filesystem::directory_iterator"
```

---

### Task 8: CometAdapter — replace QRegularExpression

**Files:**
- Modify: `src/topp/CometAdapter.cpp:38,297-302`

- [ ] **Step 1: Replace include**

Remove line 38:
```cpp
// REMOVE:
#include <QRegularExpression>
```

Add:
```cpp
#include <regex>
```

- [ ] **Step 2: Replace regex usage (lines 297-302)**

```cpp
// OLD:
    // comet_version is something like "# comet_version 2017.01 rev. 1"
    QRegularExpression comet_version_regex("(\\d{4})\\.(\\d*)rev");
    if (auto match = comet_version_regex.match(QString::fromStdString(comet_version).remove(' ')); match.hasMatch())
    {
      const int comet_year = match.captured(1).toInt();

// NEW:
    // comet_version is something like "# comet_version 2017.01 rev. 1"
    // Remove spaces for matching
    std::string version_no_spaces = comet_version;
    version_no_spaces.erase(std::remove(version_no_spaces.begin(), version_no_spaces.end(), ' '), version_no_spaces.end());
    std::regex comet_version_regex("(\\d{4})\\.(\\d*)rev");
    std::smatch match;
    if (std::regex_search(version_no_spaces, match, comet_version_regex))
    {
      const int comet_year = std::stoi(match[1].str());
```

Also add `#include <algorithm>` if not already present (for `std::remove`).

- [ ] **Step 3: Build and verify**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc)`
Expected: Compiles without Qt.

- [ ] **Step 4: Commit**

```bash
git add src/topp/CometAdapter.cpp
git commit -m "refactor(CometAdapter): replace QRegularExpression with std::regex"
```

---

### Task 9: MetaProSIP — replace all Qt usage (largest task)

**Files:**
- Modify: `src/topp/MetaProSIP.cpp:40-44,462,533-560,565,616-637,764,837-858,2038,2050-2067,2094-2114,3055-3063,3271`

This is the most complex tool change. MetaProSIP uses Qt for:
- `QProcess` — 5 instances of external R script execution
- `QFile` — copy and remove operations
- `QDir` — directory creation and absolute path
- `QFileInfo` — baseName extraction
- `QString` — function parameters and string conversions

- [ ] **Step 1: Replace Qt includes (lines 40-44)**

```cpp
// REMOVE:
#include <QtCore/QStringList>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QProcess>

// ADD:
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <filesystem>
```

- [ ] **Step 2: Create a helper function for R process execution**

MetaProSIP has 5 nearly identical QProcess blocks. Extract a shared helper at the top of the anonymous namespace or in the `RIntegration` class. Place it before `MetaProSIPReporting` (around line 460):

```cpp
/// Run an R script and return the exit status. Optionally capture merged output.
static ExternalProcess::RETURNSTATE runRScript(
    const String& r_executable,
    const String& script_path,
    const String& tmp_path,
    Size debug_level,
    String* captured_output = nullptr)
{
  std::vector<String> args = {"--vanilla"};
  if (debug_level < 1)
  {
    args.emplace_back("--quiet");
  }
  args.emplace_back("--slave");
  args.emplace_back("--file=" + script_path);

  String merged_output;
  auto capture = [&merged_output](const String& s) { merged_output += s; };

  ExternalProcess proc(capture, capture); // merge stdout+stderr into same buffer
  String error_msg;
  std::map<String, String> env = {{"R_LIBS", tmp_path}};

  auto state = proc.run(r_executable, args, "", false, error_msg, ExternalProcess::IO_MODE::READ_ONLY, env);

  if (captured_output)
  {
    *captured_output = merged_output;
  }
  return state;
}
```

- [ ] **Step 3: Change function signatures — remove QString parameters**

Replace `const QString& executable = QString("R")` with `const String& executable = "R"` in all 4 static functions:

Line 462 (`plotHeatMap`):
```cpp
// OLD:
  static void plotHeatMap(... Size debug_level = 0, const QString& executable = QString("R"))
// NEW:
  static void plotHeatMap(... Size debug_level = 0, const String& executable = "R")
```

Line 565 (`plotFilteredSpectra`):
```cpp
// OLD:
  static void plotFilteredSpectra(... Size debug_level = 0, const QString& executable = QString("R"))
// NEW:
  static void plotFilteredSpectra(... Size debug_level = 0, const String& executable = "R")
```

Line 764 (`plotScoresAndWeights`):
```cpp
// OLD:
  static void plotScoresAndWeights(... double score_plot_yaxis_min, Size debug_level = 0, const QString& executable = QString("R"))
// NEW:
  static void plotScoresAndWeights(... double score_plot_yaxis_min, Size debug_level = 0, const String& executable = "R")
```

Line 2038 (`checkRDependencies`):
```cpp
// OLD:
  static bool checkRDependencies(const String& tmp_path, StringList package_names, const QString& executable = QString("R"))
// NEW:
  static bool checkRDependencies(const String& tmp_path, StringList package_names, const String& executable = "R")
```

- [ ] **Step 4: Replace QProcess in plotHeatMap (lines 533-560)**

```cpp
// OLD (lines 533-560):
    QProcess p;
    QStringList env = QProcess::systemEnvironment();
    env << QString("R_LIBS=") + QString::fromStdString(tmp_path);
    p.setEnvironment(env);

    QStringList qparam;
    qparam << "--vanilla";
    if (debug_level < 1)
    {
      qparam << "--quiet";
    }
    qparam << "--slave" << "--file=" + QString::fromStdString(tmp_path + "/" + script_filename);
    p.start(executable, qparam);
    p.waitForFinished(-1);
    int status = p.exitCode();
    // ...
    if (status != 0)
    {
      std::cerr << "Error: Process returned with non 0 status." << std::endl;
    }
    else
    {
      QFile(QString::fromStdString(tmp_path + "/" + filename)).copy(QString::fromStdString(output_dir + "/heatmap" + file_suffix + "." + file_extension));
      if (debug_level < 1)
      {
        QFile(QString::fromStdString(tmp_path + "/" + script_filename)).remove();
        QFile(QString::fromStdString(tmp_path + "/" + filename)).remove();
      }
    }

// NEW:
    auto state = runRScript(executable, tmp_path + "/" + script_filename, tmp_path, debug_level);
    if (state != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      std::cerr << "Error: Process returned with non 0 status." << std::endl;
    }
    else
    {
      File::copy(tmp_path + "/" + filename, output_dir + "/heatmap" + file_suffix + "." + file_extension);
      if (debug_level < 1)
      {
        File::remove(tmp_path + "/" + script_filename);
        File::remove(tmp_path + "/" + filename);
      }
    }
```

- [ ] **Step 5: Replace QProcess in plotFilteredSpectra (lines 616-637)**

Same pattern as step 4. Replace the QProcess block with `runRScript` and `File::copy`/`File::remove`:

```cpp
// NEW:
      auto state = runRScript(executable, tmp_path + "/" + script_filename, tmp_path, 0); // debug_level hardcoded to 0 (quiet) in original
      if (state != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        std::cerr << "Error: Process returned with non 0 status." << std::endl;
      }
      else
      {
        File::copy(tmp_path + "/" + filename, output_dir + "/spectrum" + file_suffix + "_rt_" + String(sip_peptides[i].feature_rt) + "." + file_extension);
        if (debug_level < 1)
        {
          File::remove(tmp_path + "/" + script_filename);
          File::remove(tmp_path + "/" + filename);
        }
      }
```

Note: In the original, the QStringList always has `--quiet` (no conditional). Pass `debug_level = 0` to `runRScript` to match.

- [ ] **Step 6: Replace QProcess in plotScoresAndWeights (lines 837-858)**

Same pattern. Replace with `runRScript` and `File::copy`/`File::remove`.

- [ ] **Step 7: Replace QProcess in checkRDependencies — first block (lines 2050-2067)**

```cpp
// OLD:
      QProcess p;
      p.setProcessChannelMode(QProcess::MergedChannels);
      QStringList env = QProcess::systemEnvironment();
      env << QString("R_LIBS=") + QString::fromStdString(tmp_path);
      p.setEnvironment(env);

      QStringList checkRinPathQParam;
      checkRinPathQParam << "--vanilla" << "--quiet" << "--slave" << "--file=" + QString::fromStdString(script_filename);
      p.start(executable, checkRinPathQParam);
      p.waitForFinished(-1);

      if (p.error() == QProcess::FailedToStart || p.exitStatus() == QProcess::CrashExit || p.exitCode() != 0)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Can't execute R. ..." << std::endl;
        return false;
      }

// NEW:
      auto state = runRScript(executable, script_filename, tmp_path, 0);
      if (state != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Can't execute R. Do you have R installed? Check if the path to R is in your system path variable." << std::endl;
        return false;
      }
```

- [ ] **Step 8: Replace QProcess in checkRDependencies — second block (lines 2094-2114)**

```cpp
// OLD:
    QProcess p;
    p.setProcessChannelMode(QProcess::MergedChannels);
    QStringList env = QProcess::systemEnvironment();
    env << QString("R_LIBS=") + QString::fromStdString(tmp_path);
    p.setEnvironment(env);

    QStringList qparam;
    qparam << "--vanilla" << "--quiet" << "--slave" << "--file=" + QString::fromStdString(script_filename);
    p.start(executable, qparam);
    p.waitForFinished(-1);
    int status = p.exitCode();

    if (status != 0)
    {
      OPENMS_LOG_ERROR << "\nProblem finding all R dependencies. ..." << std::endl;
      for (TextFile::ConstIterator line_it = current_script.begin(); line_it != current_script.end(); ++line_it)
      {
        OPENMS_LOG_ERROR << *line_it  << std::endl;
      }
      QString s = p.readAllStandardOutput();
      OPENMS_LOG_ERROR << s.toStdString() << std::endl;
      return false;
    }

// NEW:
    String captured_output;
    auto state = runRScript(executable, script_filename, tmp_path, 0, &captured_output);

    if (state != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      OPENMS_LOG_ERROR << "\nProblem finding all R dependencies. Check if R and following libraries are installed:" << std::endl;
      for (TextFile::ConstIterator line_it = current_script.begin(); line_it != current_script.end(); ++line_it)
      {
        OPENMS_LOG_ERROR << *line_it  << std::endl;
      }
      OPENMS_LOG_ERROR << captured_output << std::endl;
      return false;
    }
```

- [ ] **Step 9: Replace QDir and QFileInfo in main run() function (lines 3055-3063, 3271)**

Line 3055:
```cpp
// OLD:
      QString executable = QString::fromStdString(getStringOption_("r_executable"));
      // convert path to absolute path
      QDir qc_dir(QString::fromStdString(qc_output_directory));
      qc_output_directory = qc_dir.absolutePath().toStdString();

      // trying to create qc_output_directory if not present
      if (!qc_dir.exists())
      {
        qc_dir.mkpath(QString::fromStdString(qc_output_directory));
      }

// NEW:
      String executable = getStringOption_("r_executable");
      // convert path to absolute path
      qc_output_directory = File::absolutePath(qc_output_directory);

      // trying to create qc_output_directory if not present
      if (!File::exists(qc_output_directory))
      {
        std::filesystem::create_directories(std::string(qc_output_directory));
      }
```

Line 3271:
```cpp
// OLD:
    String file_suffix = "_" + String(QFileInfo(QString::fromStdString(in_mzml)).baseName().toStdString()) + "_" + String::random(4);

// NEW:
    String file_suffix = "_" + FileHandler::stripExtension(File::basename(in_mzml)) + "_" + String::random(4);
```

- [ ] **Step 10: Build and verify**

Run: `cmake --build OpenMS-build --target MetaProSIP -j$(nproc)`
Expected: Compiles without Qt. No Qt headers remain.

- [ ] **Step 11: Verify no remaining Qt references**

Run: `grep -n "Q[A-Z]" src/topp/MetaProSIP.cpp` (should return nothing)

- [ ] **Step 12: Commit**

```bash
git add src/topp/MetaProSIP.cpp
git commit -m "refactor(MetaProSIP): replace QProcess/QFile/QDir with ExternalProcess and File utilities"
```

---

### Task 10: Phase 1 CMake cleanup — remove Qt6::Core linking

**Files:**
- Modify: `src/topp/CMakeLists.txt:73-85`

- [ ] **Step 1: Remove the Qt6::Core linking block**

In `src/topp/CMakeLists.txt`, remove lines 73-85:

```cmake
# REMOVE this entire block:
# Qt6::Core for tools that directly use Qt headers (optional when Qt not available)
if(TARGET Qt6::Core)
  target_link_libraries(AssayGeneratorMetaboSirius Qt6::Core)
  target_link_libraries(CometAdapter Qt6::Core)
  target_link_libraries(IDRipper Qt6::Core)
  target_link_libraries(MetaProSIP Qt6::Core)
  target_link_libraries(MzMLSplitter Qt6::Core)
  target_link_libraries(OpenSwathFileSplitter Qt6::Core)
  target_link_libraries(OpenSwathWorkflow Qt6::Core)
  target_link_libraries(QCEmbedder Qt6::Core)
  target_link_libraries(QCExtractor Qt6::Core)
  target_link_libraries(QCShrinker Qt6::Core)
endif()
```

- [ ] **Step 2: Full build test**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: All TOPP tools compile without Qt6::Core.

- [ ] **Step 3: Verify no Qt includes remain in src/topp/**

Run: `grep -rn "#include.*<Q" src/topp/ --include="*.cpp"` — should only match the 3 GUI tools (ExecutePipeline, INIUpdater, ImageCreator) via their OpenMS_GUI includes.

Actually, also check for `QtCore`, `QtGui` patterns:
Run: `grep -rn "Qt[A-Z]" src/topp/ --include="*.cpp"` — should only match the 3 GUI tools.

- [ ] **Step 4: Commit**

```bash
git add src/topp/CMakeLists.txt
git commit -m "build(topp): remove Qt6::Core linking — all non-GUI TOPP tools are Qt-free"
```

---

### Task 11: Move GUI tools to openms_gui (Phase 2)

**Files:**
- Move: `src/topp/ExecutePipeline.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/`
- Move: `src/topp/INIUpdater.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/`
- Move: `src/topp/ImageCreator.cpp` → `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/`
- Modify: `src/topp/executables.cmake:161-168`
- Modify: `src/topp/CMakeLists.txt:53-64`
- Modify: `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/executables.cmake:5-11`
- Modify: `src/openms/source/APPLICATIONS/ToolHandler.cpp:196-207`

- [ ] **Step 1: Move source files**

```bash
git mv src/topp/ExecutePipeline.cpp src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/ExecutePipeline.cpp
git mv src/topp/INIUpdater.cpp src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/INIUpdater.cpp
git mv src/topp/ImageCreator.cpp src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/ImageCreator.cpp
```

- [ ] **Step 2: Remove TOPP_executables_with_GUIlib from executables.cmake**

In `src/topp/executables.cmake`, remove lines 161-168 entirely:

```cmake
# REMOVE:
## all targets requiring OpenMS_GUI
set(TOPP_executables_with_GUIlib
ExecutePipeline
# util category
ImageCreator
INIUpdater
)
```

Also update the `foreach` on line 172 to remove the reference:
```cmake
# OLD:
foreach(i ${TOPP_executables} ${TOPP_executables_with_GUIlib})

# NEW:
foreach(i ${TOPP_executables})
```

- [ ] **Step 3: Remove WITH_GUI block from topp/CMakeLists.txt**

Remove lines 53-64:
```cmake
# REMOVE:
# some regular TOPP tools need the GUI lib, only build them when WITH_GUI is enabled
if(WITH_GUI)
  foreach(i ${TOPP_executables_with_GUIlib})
    add_executable(${i} ${i}.cpp)
    # we also want to install each topp tool
    openms_add_executable_compiler_flags(${i})
    install_tool(${i})
    target_link_libraries(${i} OpenMS OpenMS_GUI)
  endforeach(i)

  set(TOPP_executables ${TOPP_executables} ${TOPP_executables_with_GUIlib})
endif()
```

- [ ] **Step 4: Add tools to GUI executables.cmake**

In `src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/executables.cmake`, add the 3 tools to the `GUI_executables` list (alphabetically):

```cmake
set(GUI_executables
ExecutePipeline
FLASHDeconvWizard
ImageCreator
INIFileEditor
INIUpdater
SwathWizard
TOPPAS
TOPPView
)
```

- [ ] **Step 5: Ensure GUI tools are in TOPP_TOOLS for CWL/INI generation**

The 3 moved tools must still appear in the `TOPP_TOOLS` cache variable so that CWL generation and `-write_ini` tests work. Add to the bottom of `src/openms_gui/CMakeLists.txt` (after the GUI executable build loop, around line 91):

```cmake
# Add GUI-based TOPP tools to the TOPP_TOOLS list for CWL/INI generation
set(GUI_TOPP_TOOLS ExecutePipeline ImageCreator INIUpdater)
set(TOPP_TOOLS ${TOPP_TOOLS} ${GUI_TOPP_TOOLS}
    CACHE INTERNAL "OpenMS' TOPP tools" FORCE)
```

Also add them to the `TOPP` collection target dependency. In `src/openms_gui/CMakeLists.txt`, after the executable build loop:
```cmake
# Make the TOPP collection target depend on GUI TOPP tools
if(TARGET TOPP)
  add_dependencies(TOPP ${GUI_TOPP_TOOLS})
endif()
```

- [ ] **Step 6: Verify linking in GUI CMakeLists**

The existing GUI build loop (line 90) does `target_link_libraries(${i} OpenMS_GUI)`. Since `OpenMS_GUI` links `OpenMS` as a public dependency, the 3 moved tools get `OpenMS` transitively. Verify this is sufficient — ExecutePipeline, INIUpdater, and ImageCreator previously linked both `OpenMS` and `OpenMS_GUI` explicitly. If transitive linking isn't enough (compilation succeeds but linking fails), add:

```cmake
# After the GUI executable loop (line 91), add explicit OpenMS link for TOPP-origin tools:
foreach(i ExecutePipeline ImageCreator INIUpdater)
  target_link_libraries(${i} OpenMS)
endforeach()
```

- [ ] **Step 7: Build with WITH_GUI=ON**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: All tools compile. ExecutePipeline, INIUpdater, ImageCreator built from new location.

- [ ] **Step 8: Build with WITH_GUI=OFF (if feasible to test)**

Run: `cmake -DWITH_GUI=OFF OpenMS-build && cmake --build OpenMS-build -j$(nproc)`
Expected: All regular TOPP tools compile. No Qt dependency. GUI tools skipped.

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "refactor: move ExecutePipeline, INIUpdater, ImageCreator to openms_gui

These 3 TOPP tools depend on OpenMS_GUI (TOPPASScene, QImage/QPainter).
Moving them to src/openms_gui/ makes src/topp/ fully Qt-free.
They remain in TOPP_TOOLS for CWL/INI generation."
```

---

### Task 12: CMake cleanup — move Qt find inside WITH_GUI, guard KNIME

**Files:**
- Modify: `cmake/cmake_findExternalLibs.cmake:342-350`
- Modify: `cmake/knime_package_support.cmake:68-71`

- [ ] **Step 1: Move Qt6::Core find_package inside WITH_GUI block**

In `cmake/cmake_findExternalLibs.cmake`, replace lines 342-350:

```cmake
# OLD:
# Qt6::Core is needed by TOPP tools and GUI, but NOT by libOpenMS itself.
# When building without GUI (e.g., pyOpenMS wheels), Qt is optional.
find_package(Qt6 ${QT_MIN_VERSION} COMPONENTS Core QUIET)

IF (Qt6Core_FOUND)
  message(STATUS "Found Qt ${Qt6Core_VERSION}")
ELSE()
  message(STATUS "Qt6Core not found — TOPP tools and GUI will not be available")
ENDIF()
```

Move this into the existing `if (WITH_GUI)` block (line 363), placing it before the additional component search. The new block at line 363 becomes:

```cmake
if (WITH_GUI)
  # Qt6::Core is needed by GUI library and GUI-based TOPP tools
  find_package(Qt6 ${QT_MIN_VERSION} COMPONENTS Core QUIET)

  IF (Qt6Core_FOUND)
    message(STATUS "Found Qt ${Qt6Core_VERSION}")
  ELSE()
    message(FATAL_ERROR "Qt6Core not found — required when WITH_GUI=ON. Use -DWITH_GUI=OFF to build without GUI.")
  ENDIF()

  # --------------------------------------------------------------------------
  # Find additional Qt libs
  # ... (rest of existing WITH_GUI block unchanged)
```

Note: Changed from `message(STATUS)` to `message(FATAL_ERROR)` when Qt is not found AND WITH_GUI=ON, since GUI build now requires it.

- [ ] **Step 2: Guard KNIME Qt usage**

In `cmake/knime_package_support.cmake`, guard lines 68-71:

```cmake
# OLD:
# Qt6::Core is already found by cmake_findExternalLibs.cmake
get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_LIBS" OUTPUT_VARIABLE QT_INSTALL_LIBS)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_BINS" OUTPUT_VARIABLE QT_INSTALL_BINS)

# NEW:
if(TARGET Qt6::qmake)
  get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
  exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_LIBS" OUTPUT_VARIABLE QT_INSTALL_LIBS)
  exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_BINS" OUTPUT_VARIABLE QT_INSTALL_BINS)
endif()
```

- [ ] **Step 3: Build and verify**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Build succeeds with Qt present.

- [ ] **Step 4: Commit**

```bash
git add cmake/cmake_findExternalLibs.cmake cmake/knime_package_support.cmake
git commit -m "build: move Qt6 find_package inside WITH_GUI block

Qt is no longer needed by any non-GUI component. Guard KNIME
packaging against missing Qt targets."
```

---

### Task 13: Final verification

- [ ] **Step 1: Verify no Qt references in src/topp/**

```bash
grep -rn "#include.*<Q" src/topp/ --include="*.cpp"
grep -rn "Qt[A-Z]" src/topp/ --include="*.cpp"
grep -rn "QString\|QFile\|QDir\|QProcess\|QRegular" src/topp/ --include="*.cpp"
```

All three should return zero matches.

- [ ] **Step 2: Run full test suite**

Run: `ctest --test-dir OpenMS-build -j$(nproc)`
Expected: All tests pass (same as before the refactor).

- [ ] **Step 3: Verify CWL and INI generation still works**

Run: `cmake --build OpenMS-build --target generate_cwl_files` (if available)
Check that ExecutePipeline, INIUpdater, ImageCreator still appear in generated artifacts.

- [ ] **Step 4: Verify Resampler tests run without WITH_GUI gate**

Run: `ctest --test-dir OpenMS-build -R "TOPP_Resampler" -V`
Expected: Both Resampler tests pass.
