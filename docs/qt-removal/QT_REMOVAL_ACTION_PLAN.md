# Qt Removal from OpenMS Core - Action Plan

## Issue Reference
GitHub Issue: [#8200 - Remove all Qt dependencies from the OpenMS and OpenSwath libraries](https://github.com/OpenMS/OpenMS/issues/8200)

## Objective
Remove all Qt dependencies from the OpenMS core library and OpenSwath, while keeping Qt support in the GUI library (`openms_gui`) and GUI executables.

## Current Status (Based on Audit)

### Statistics
- **60 files** need Qt removal
- **134 QString occurrences**
- **61 QDir occurrences**
- **48 QStringList occurrences**
- **30 QDateTime occurrences**
- **29 QProcess occurrences**
- **26 QFileInfo occurrences**
- **25 QByteArray occurrences**
- **18 QFile occurrences**
- **15 QUrl occurrences**
- **7 QNetworkAccessManager occurrences**
- **7 QFileSystemWatcher occurrences**

### Key Files to Modify

#### Core Data Structures (High Priority - Wide Impact)
1. `src/openms/include/OpenMS/DATASTRUCTURES/String.h` & `.cpp`
2. `src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h` & `.cpp`
3. `src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h`
4. `src/openms/include/OpenMS/DATASTRUCTURES/StringListUtils.h` & `.cpp`
5. `src/openms/include/OpenMS/DATASTRUCTURES/DateTime.cpp`

#### System Utilities
6. `src/openms/include/OpenMS/SYSTEM/File.h` & `.cpp` (QFile, QDir, QFileInfo)
7. `src/openms/include/OpenMS/SYSTEM/ExternalProcess.h` & `.cpp` (QProcess)
8. `src/openms/include/OpenMS/SYSTEM/NetworkGetRequest.h` & `.cpp` (QNetworkAccessManager)
9. `src/openms/include/OpenMS/SYSTEM/FileWatcher.h` & `.cpp` (QFileSystemWatcher)
10. `src/openms/include/OpenMS/SYSTEM/RWrapper.h` & `.cpp`
11. `src/openms/include/OpenMS/SYSTEM/PythonInfo.cpp`
12. `src/openms/include/OpenMS/SYSTEM/JavaInfo.cpp`
13. `src/openms/include/OpenMS/SYSTEM/UpdateCheck.cpp`

#### Application Framework
14. `src/openms/include/OpenMS/APPLICATIONS/TOPPBase.h` & `.cpp`
15. `src/openms/include/OpenMS/APPLICATIONS/ToolHandler.h` & `.cpp`

#### Format Handlers
16. `src/openms/source/FORMAT/MascotRemoteQuery.h` & `.cpp`
17. Other format files using QString/QFile

#### CMake
18. `src/openms/CMakeLists.txt` - Remove Qt6::Core and Qt6::Network

## Implementation Strategy

### Phase 1: Core Data Structure Refactoring (Week 1-2)

#### 1.1 String Class Refactoring
**Files**: `String.h`, `String.cpp`, `StringUtils.h`, `StringListUtils.h/cpp`

**Changes**:
- Remove `String::toQString()` method
- Remove `QString` constructor from `String`
- Remove `StringListUtils::fromQStringList()`
- Add helper utilities for boundary conversion (GUI still needs Qt)

**Testing**: Run all unit tests in `String_test.cpp`

#### 1.2 DataValue Refactoring
**Files**: `DataValue.h`, `DataValue.cpp`

**Changes**:
- Remove `DataValue::toQString()` method
- Remove `DataValue(const QString&)` constructor
- Use `std::string` internally

**Testing**: Run `DataValue_test.cpp`

#### 1.3 DateTime Refactoring
**Files**: `DateTime.cpp`

**Changes**:
- Replace `QDateTime` with `std::chrono` + custom formatting
- Consider using `boost::date_time` if complex formatting needed

**Testing**: Run `DateTime_test.cpp`

### Phase 2: File System Operations (Week 2-3)

#### 2.1 File Class Refactoring
**Files**: `File.h`, `File.cpp`

**Changes**:
- Replace `QFile`, `QFileInfo`, `QDir` with `std::filesystem` (C++17)
- Implement:
  - `exists()` using `std::filesystem::exists()`
  - `readable()` using `std::filesystem::perms`
  - `writable()` using `std::filesystem::perms`
  - `isDirectory()` using `std::filesystem::is_directory()`
  - File size using `std::filesystem::file_size()`

**Example**:
```cpp
// Before:
QFileInfo fi(file.toQString());
return fi.exists();

// After:
std::filesystem::path p(file);
return std::filesystem::exists(p);
```

**Testing**: Run `File_test.cpp`

### Phase 3: Process Management (Week 3-4)

#### 3.1 ExternalProcess Refactoring
**Files**: `ExternalProcess.h`, `ExternalProcess.cpp`

**Changes**:
- Replace `QProcess` with:
  - Option 1: `boost::process` (recommended if already using Boost)
  - Option 2: Custom implementation using `fork()` + `exec()` on Linux/macOS, `CreateProcess()` on Windows
  - Option 3: `std::system()` for simple cases (not recommended)

**Testing**: Run `ExternalProcess_test.cpp`

#### 3.2 RWrapper, PythonInfo, JavaInfo
**Files**: Multiple `*Info.cpp` files

**Changes**:
- Use new `ExternalProcess` API
- Replace `QString`/`QStringList` with `std::string`/`std::vector<std::string>`

### Phase 4: Networking (Week 4-5)

#### 4.1 NetworkGetRequest Refactoring
**Files**: `NetworkGetRequest.h`, `NetworkGetRequest.cpp`

**Note**: PR #8201 is already addressing this with `cpp-httplib`. Check its status first.

**Changes**:
- Replace `QNetworkAccessManager` with `cpp-httplib` or `boost::asio`
- Update `UpdateCheck.cpp` to use new API

**Testing**: Run network-related tests

### Phase 5: File Watching (Week 5)

#### 5.1 FileWatcher Refactoring
**Files**: `FileWatcher.h`, `FileWatcher.cpp`

**Changes**:
- Replace `QFileSystemWatcher` with:
  - Linux: `inotify` API
  - Windows: `ReadDirectoryChangesW`
  - macOS: `FSEvents`
  - Or use cross-platform library like `efsw` or `boost::asio::file_descriptor`

**Testing**: Run `FileWatcher_test.cpp`

### Phase 6: Application Framework (Week 6)

#### 6.1 TOPPBase and ToolHandler
**Files**: `TOPPBase.h/cpp`, `ToolHandler.h/cpp`

**Changes**:
- Replace `QString` with `std::string`
- Update parameter handling

**Testing**: Run TOPP tool tests

### Phase 7: CMake Configuration (Week 6)

#### 7.1 Remove Qt from CMakeLists.txt
**File**: `src/openms/CMakeLists.txt`

**Changes**:
```cmake
# Before:
set(OPENMS_DEP_LIBRARIES Evergreen LibSVM::LibSVM XercesC::XercesC Qt6::Core Qt6::Network)

# After:
set(OPENMS_DEP_LIBRARIES Evergreen LibSVM::LibSVM XercesC::XercesC)
```

**Testing**: Build entire project and verify `openms_gui` still works

### Phase 8: GUI Boundary Layer (Week 7)

#### 8.1 Create Conversion Utilities
Since `openms_gui` still needs Qt, create boundary conversion utilities:

**New file**: `src/openms/include/OpenMS/DATASTRUCTURES/QtConversion.h`

```cpp
#ifdef WITH_GUI
#include <QString>
#include <QStringList>

namespace OpenMS {
  inline QString toQt(const std::string& s) { return QString::fromStdString(s); }
  inline std::string fromQt(const QString& s) { return s.toStdString(); }
  inline QStringList toQt(const std::vector<std::string>& v) {
    QStringList result;
    for (const auto& s : v) result << QString::fromStdString(s);
    return result;
  }
}
#endif
```

### Phase 9: Testing & Validation (Week 7-8)

1. **Build Tests**: Ensure project builds without Qt for core
2. **Unit Tests**: Run full test suite
3. **pyOpenMS**: Verify Python bindings work
4. **GUI Tests**: Ensure `openms_gui` still functions
5. **Performance**: Benchmark to ensure no regressions

## Git Workflow

### 1. Create Feature Branch
```bash
cd /home/prachi/openMS/OpenMS
git checkout develop
git pull upstream develop
git checkout -b feature/remove-qt-from-core
```

### 2. Commit Strategy
Make small, focused commits:
```bash
git commit -m "Refactor: Remove QString from String class"
git commit -m "Refactor: Replace QFile with std::filesystem in File class"
git commit -m "Refactor: Replace QProcess in ExternalProcess"
git commit -m "Build: Remove Qt dependencies from OpenMS CMakeLists.txt"
```

### 3. Testing Between Commits
```bash
cmake -B build -S . -DWITH_GUI=OFF -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target OpenMS -j4
ctest --test-dir build -R OpenMS
```

### 4. Create Pull Request
```bash
git push origin feature/remove-qt-from-core
```

Then create PR on GitHub with:
- Title: "Remove all Qt dependencies from OpenMS and OpenSwath libraries"
- Body: "Fixes #8200"
- Detailed description of changes

## Comment to Post on GitHub Issue

```markdown
Hi @jpfeuffer and @timosachsenberg,

I'd like to work on this issue! I've done an audit of the codebase and created an implementation plan.

**Current Status:**
- 60 files need Qt removal
- Main Qt usage: QString (134), QDir (61), QStringList (48), QDateTime (30), QProcess (29)
- OpenSwathAlgo is already Qt-free ✓

**Proposed Approach:**
1. Phase 1: Core data structures (String, DataValue, DateTime)
2. Phase 2: File system operations (File.cpp using std::filesystem)
3. Phase 3: Process management (ExternalProcess using boost::process)
4. Phase 4: Networking (coordinating with PR #8201)
5. Phase 5-7: FileWatcher, TOPP tools, CMake updates
6. Phase 8: GUI boundary layer with conversion helpers

I'll submit small, focused PRs for each phase with comprehensive tests. I have a Linux dev environment set up and have already built OpenMS successfully.

Would this approach work for the maintainers? Any concerns or suggestions?
```

## Dependencies & Libraries to Consider

- **C++17/20 Standard Library**: `std::filesystem`, `std::chrono`, `std::string`
- **Boost** (already used in OpenMS):
  - `boost::process` for external processes
  - `boost::date_time` if chrono is insufficient
  - `boost::filesystem` as fallback
- **cpp-httplib**: For networking (check PR #8201)
- **Platform APIs**: For FileWatcher (inotify on Linux)

## Success Criteria

- [ ] All Qt includes removed from `src/openms/` (except extern)
- [ ] All Qt includes removed from `src/openswathalgo/`
- [ ] `src/openms/CMakeLists.txt` no longer links Qt
- [ ] All unit tests pass
- [ ] pyOpenMS builds successfully
- [ ] GUI library and tools still work with Qt
- [ ] Build time improved (target: < 2 hours)
- [ ] No functionality regression

## Resources

- Issue: https://github.com/OpenMS/OpenMS/issues/8200
- Related PR (networking): https://github.com/OpenMS/OpenMS/pull/8201
- C++ filesystem: https://en.cppreference.com/w/cpp/filesystem
- Boost.Process: https://www.boost.org/doc/libs/1_83_0/doc/html/process.html
