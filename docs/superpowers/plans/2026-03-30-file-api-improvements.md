# File API Improvements — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `stemName()`, `extension()`, and `listDirectories()` to `OpenMS::File`, then replace awkward patterns in TOPP tools.

**Architecture:** Add 3 new static methods to `File` in the core library. `stemName()` delegates to `FileHandler::stripExtension(File::basename())`. `extension()` extracts compound-aware extensions. `listDirectories()` wraps `std::filesystem::directory_iterator` with directory filtering. Then replace all call sites in `src/topp/` and fix missed `File::fileSize()` / `File::makeDir()` opportunities.

**Tech Stack:** C++20, OpenMS::File, OpenMS::FileHandler, std::filesystem

---

## File Map

- Modify: `src/openms/include/OpenMS/SYSTEM/File.h` — add 3 method declarations
- Modify: `src/openms/source/SYSTEM/File.cpp` — add 3 method implementations
- Modify: `src/tests/class_tests/openms/source/File_test.cpp` — add tests for all 3 methods
- Modify: 11 TOPP tool `.cpp` files — replace patterns with new APIs

---

### Task 1: Add `File::stemName()` with tests

**Files:**
- Modify: `src/openms/include/OpenMS/SYSTEM/File.h`
- Modify: `src/openms/source/SYSTEM/File.cpp`
- Modify: `src/tests/class_tests/openms/source/File_test.cpp`

- [ ] **Step 1: Add declaration to File.h**

After the `basename()` declaration (line 138), add:

```cpp
    /// Returns the basename of the file without any known file extension.
    /// Delegates to FileHandler::stripExtension(File::basename(file)).
    /// E.g., "/path/sample.mzML.gz" returns "sample", "/path/data.featureXML" returns "data".
    /// Unknown extensions are stripped at the last dot: "/path/file.txt" returns "file".
    /// Directories with dots in the path are handled correctly: "/my.dir/file" returns "file".
    static String stemName(const String& file);
```

- [ ] **Step 2: Add implementation to File.cpp**

Add the include at the top of File.cpp (check if already present):
```cpp
#include <OpenMS/FORMAT/FileHandler.h>
```

Add the implementation after `basename()` (after line 333):

```cpp
String File::stemName(const String& file)
{
  return FileHandler::stripExtension(basename(file));
}
```

- [ ] **Step 3: Add tests to File_test.cpp**

Find the `basename()` test section (around line 116-120). Add a new section after it:

```cpp
START_SECTION((static String stemName(const String &file)))
  // basic: strips known extension from full path
  TEST_EQUAL(File::stemName("/path/to/sample.mzML"), "sample");
  // compound extension: .mzML.gz is a known compound extension
  TEST_EQUAL(File::stemName("/path/to/sample.mzML.gz"), "sample");
  // unknown extension: strips last dot segment
  TEST_EQUAL(File::stemName("/path/to/file.txt"), "file");
  // unknown compound: only strips known part
  TEST_EQUAL(File::stemName("/path/to/file.txt.tgz"), "file.txt");
  // no extension
  TEST_EQUAL(File::stemName("/path/to/file"), "file");
  // filename only (no path)
  TEST_EQUAL(File::stemName("experiment.featureXML"), "experiment");
  // empty string
  TEST_EQUAL(File::stemName(""), "");
  // dotted directory, no extension on file
  TEST_EQUAL(File::stemName("/home.with.dot/filename"), "filename");
  // Windows path
  TEST_EQUAL(File::stemName("c:\\data\\sample.idXML"), "sample");
  // extension-only name
  TEST_EQUAL(File::stemName(".mzML"), "");
END_SECTION
```

- [ ] **Step 4: Build and run tests**

```bash
cmake --build OpenMS-build --target File_test -j$(nproc)
OpenMS-build/bin/File_test
```

Expected: All tests pass, including the new `stemName` section.

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/SYSTEM/File.h src/openms/source/SYSTEM/File.cpp src/tests/class_tests/openms/source/File_test.cpp
git commit -m "feat(File): add stemName() — basename without known file extension"
```

---

### Task 2: Add `File::extension()` with tests

**Files:**
- Modify: `src/openms/include/OpenMS/SYSTEM/File.h`
- Modify: `src/openms/source/SYSTEM/File.cpp`
- Modify: `src/tests/class_tests/openms/source/File_test.cpp`

- [ ] **Step 1: Add declaration to File.h**

After the `stemName()` declaration, add:

```cpp
    /// Returns the file extension including the leading dot.
    /// Recognizes compound OpenMS extensions like ".mzML.gz".
    /// E.g., "/path/sample.mzML.gz" returns ".mzML.gz", "/path/file.txt" returns ".txt".
    /// Returns empty string if there is no extension: "/path/file" returns "".
    static String extension(const String& file);
```

- [ ] **Step 2: Add implementation to File.cpp**

After `stemName()`:

```cpp
String File::extension(const String& file)
{
  String base = basename(file);
  String stem = FileHandler::stripExtension(base);
  if (stem.size() >= base.size())
  {
    return ""; // no extension (stripExtension returned the same or longer string)
  }
  return base.substr(stem.size()); // everything after the stem, including leading '.'
}
```

- [ ] **Step 3: Add tests to File_test.cpp**

After the `stemName` test section:

```cpp
START_SECTION((static String extension(const String &file)))
  // known extension
  TEST_EQUAL(File::extension("/path/to/sample.mzML"), ".mzML");
  // compound extension
  TEST_EQUAL(File::extension("/path/to/sample.mzML.gz"), ".mzML.gz");
  // unknown extension
  TEST_EQUAL(File::extension("/path/to/file.txt"), ".txt");
  // no extension
  TEST_EQUAL(File::extension("/path/to/file"), "");
  // filename only
  TEST_EQUAL(File::extension("experiment.featureXML"), ".featureXML");
  // empty string
  TEST_EQUAL(File::extension(""), "");
  // dotted directory, no extension on file
  TEST_EQUAL(File::extension("/home.with.dot/filename"), "");
  // Windows path
  TEST_EQUAL(File::extension("c:\\data\\sample.idXML"), ".idXML");
  // extension-only name
  TEST_EQUAL(File::extension(".mzML"), ".mzML");
END_SECTION
```

- [ ] **Step 4: Build and run tests**

```bash
cmake --build OpenMS-build --target File_test -j$(nproc)
OpenMS-build/bin/File_test
```

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/SYSTEM/File.h src/openms/source/SYSTEM/File.cpp src/tests/class_tests/openms/source/File_test.cpp
git commit -m "feat(File): add extension() — compound-aware file extension extraction"
```

---

### Task 3: Add `File::listDirectories()` with tests

**Files:**
- Modify: `src/openms/include/OpenMS/SYSTEM/File.h`
- Modify: `src/openms/source/SYSTEM/File.cpp`
- Modify: `src/tests/class_tests/openms/source/File_test.cpp`

- [ ] **Step 1: Add declaration to File.h**

After the `extension()` declaration, add:

```cpp
    /// Returns a sorted list of subdirectory paths (non-recursive) in the given directory.
    /// Returns full absolute-style paths using '/' separators.
    /// If the path does not exist or is not a directory, returns an empty list (no throw).
    static StringList listDirectories(const String& dir);
```

- [ ] **Step 2: Add implementation to File.cpp**

After `extension()`:

```cpp
StringList File::listDirectories(const String& dir)
{
  StringList result;
  auto dir_path = to_path(dir);
  std::error_code ec;
  if (!fs::is_directory(dir_path, ec)) return result;

  for (const auto& entry : fs::directory_iterator(dir_path, ec))
  {
    if (entry.is_directory())
    {
      result.push_back(entry.path().generic_string());
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}
```

This follows the same pattern as `fileList()` (line 433 of File.cpp): uses `std::error_code` overloads, sorts results, uses `generic_string()` for consistent `/` separators.

- [ ] **Step 3: Add tests to File_test.cpp**

After the `extension` test section:

```cpp
START_SECTION((static StringList listDirectories(const String &dir)))
  // create temp structure with subdirectories
  File::TempDir tdir;
  String base = tdir.getPath();
  File::makeDir(base + "/subA");
  File::makeDir(base + "/subB");
  // also create a file (should NOT appear in results)
  {
    std::ofstream f(std::string(base + "/afile.txt"));
    f << "test";
  }

  StringList dirs = File::listDirectories(base);
  TEST_EQUAL(dirs.size(), 2);
  // results are sorted
  TEST_TRUE(dirs[0].hasSuffix("subA"));
  TEST_TRUE(dirs[1].hasSuffix("subB"));

  // non-existent directory returns empty list
  StringList empty = File::listDirectories("/nonexistent_path_xyz");
  TEST_EQUAL(empty.size(), 0);

  // empty string
  StringList from_empty = File::listDirectories("");
  // behavior: empty string → not a directory → empty list
  // (unless cwd has subdirs, but we just test it doesn't crash)
END_SECTION
```

- [ ] **Step 4: Build and run tests**

```bash
cmake --build OpenMS-build --target File_test -j$(nproc)
OpenMS-build/bin/File_test
```

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/openms/include/OpenMS/SYSTEM/File.h src/openms/source/SYSTEM/File.cpp src/tests/class_tests/openms/source/File_test.cpp
git commit -m "feat(File): add listDirectories() — sorted list of subdirectories"
```

---

### Task 4: Replace patterns in TOPP tools

**Files to modify (17 `stemName` replacements + 2 fixes):**
- `src/topp/QCEmbedder.cpp:135`
- `src/topp/QCExtractor.cpp:109`
- `src/topp/QCShrinker.cpp:107`
- `src/topp/MetaProSIP.cpp:3244`
- `src/topp/OpenSwathFileSplitter.cpp:98`
- `src/topp/MascotAdapterOnline.cpp:266`
- `src/topp/ConsensusID.cpp:687,754`
- `src/topp/ProteinQuantifier.cpp:386,571,634,773`
- `src/topp/ProteomicsLFQ.cpp:851,852,1129,1130`
- `src/topp/SpectraSTSearchAdapter.cpp:293`
- `src/topp/AssayGeneratorMetaboSirius.cpp:194` (listDirectories)
- `src/topp/MzMLSplitter.cpp:99` (File::fileSize)
- `src/topp/MetaProSIP.cpp:3036` (File::makeDir)

- [ ] **Step 1: Replace `FileHandler::stripExtension(File::basename(...))` with `File::stemName(...)` in all 17 TOPP tool call sites**

For each file, the replacement is mechanical:

```cpp
// OLD:
FileHandler::stripExtension(File::basename(some_path))

// NEW:
File::stemName(some_path)
```

In files where `FileHandler.h` was only included for `stripExtension`, remove the include. Check each file: if `FileHandler` is still used for other calls (e.g., `loadExperiment`, `storeExperiment`, `swapExtension`, `getTypeByFileName`), keep the include.

Files where `FileHandler.h` can likely be removed (only used for stripExtension):
- `QCEmbedder.cpp` — check, may still use FileHandler for loading
- `QCExtractor.cpp` — check
- `QCShrinker.cpp` — check
- `OpenSwathFileSplitter.cpp` — check

Files where `FileHandler.h` is definitely still needed:
- `MetaProSIP.cpp` — uses FileHandler for loading
- `ConsensusID.cpp` — uses FileHandler for loading
- `ProteinQuantifier.cpp` — uses FileHandler
- `ProteomicsLFQ.cpp` — uses FileHandler for loading
- `MascotAdapterOnline.cpp` — uses FileHandler
- `SpectraSTSearchAdapter.cpp` — uses FileHandler

- [ ] **Step 2: Replace `std::filesystem::directory_iterator` in AssayGeneratorMetaboSirius**

In `src/topp/AssayGeneratorMetaboSirius.cpp`, replace the directory iteration block (around line 192-200):

```cpp
// OLD:
    std::vector<String> subdirs;
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(std::string(sirius_project_directory), ec))
    {
      if (entry.is_directory())
      {
        subdirs.emplace_back(entry.path().string());
      }
    }

// NEW:
    std::vector<String> subdirs_list = File::listDirectories(sirius_project_directory);
    std::vector<String> subdirs(subdirs_list.begin(), subdirs_list.end());
```

Note: `File::listDirectories` returns `StringList` (which is `std::vector<String>`), so this can be simplified to:

```cpp
    std::vector<String> subdirs = File::listDirectories(sirius_project_directory);
```

Also remove `#include <filesystem>` if no longer needed (check for other `std::filesystem` usage in the file). Add `#include <OpenMS/SYSTEM/File.h>` if not present.

- [ ] **Step 3: Replace `std::filesystem::file_size` in MzMLSplitter**

In `src/topp/MzMLSplitter.cpp` (line 99):

```cpp
// OLD:
      float total_size = static_cast<float>(std::filesystem::file_size(std::string(in)));

// NEW:
      float total_size = static_cast<float>(File::fileSize(in));
```

Remove `#include <filesystem>`. Add `#include <OpenMS/SYSTEM/File.h>` if not present.

- [ ] **Step 4: Replace `std::filesystem::create_directories` in MetaProSIP**

In `src/topp/MetaProSIP.cpp` (around line 3036):

```cpp
// OLD:
        std::filesystem::create_directories(std::string(qc_output_directory));

// NEW:
        File::makeDir(qc_output_directory);
```

Check if `#include <filesystem>` can be removed (search for other `std::filesystem` usage in the file).

- [ ] **Step 5: Build all modified tools**

```bash
cmake --build OpenMS-build -j$(nproc)
```

Expected: Full build succeeds.

- [ ] **Step 6: Run tests for modified tools**

```bash
ctest --test-dir OpenMS-build -R "TOPP_(QCEmbedder|QCExtractor|QCShrinker|MetaProSIP|OpenSwathFileSplitter|MascotAdapterOnline|ConsensusID|ProteinQuantifier|ProteomicsLFQ|SpectraSTSearchAdapter|AssayGeneratorMetabo|MzMLSplitter)" --output-on-failure
```

Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/topp/*.cpp
git commit -m "refactor(topp): use File::stemName(), File::listDirectories(), File::fileSize(), File::makeDir()

Replace FileHandler::stripExtension(File::basename()) with File::stemName() in 17 call sites.
Replace std::filesystem::directory_iterator with File::listDirectories().
Replace std::filesystem::file_size with File::fileSize().
Replace std::filesystem::create_directories with File::makeDir()."
```
