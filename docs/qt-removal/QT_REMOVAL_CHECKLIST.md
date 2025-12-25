# Qt Removal Checklist - Issue #8200

## Pre-Work
- [x] Fork OpenMS repository
- [x] Clone to local machine
- [x] Create audit script
- [x] Run audit and identify all Qt usage
- [ ] Build OpenMS successfully (baseline)
- [ ] Create feature branch: `feature/remove-qt-from-core`
- [ ] Post plan on GitHub issue

## Phase 1: Core Data Structures

### String & StringUtils
- [ ] Remove `String::toQString()` method
- [ ] Remove `String(const QString&)` constructor
- [ ] Update `StringUtils` to remove Qt methods
- [ ] Update `StringListUtils::fromQStringList()`
- [ ] Run `String_test.cpp`
- [ ] Commit: "Refactor: Remove QString from String class"

### DataValue
- [ ] Remove `DataValue::toQString()` method
- [ ] Remove `DataValue(const QString&)` constructor
- [ ] Remove `operator=(const QString&)`
- [ ] Run `DataValue_test.cpp`
- [ ] Commit: "Refactor: Remove QString from DataValue"

### DateTime
- [ ] Replace QDateTime with std::chrono
- [ ] Implement custom date formatting
- [ ] Update `toString()` and `fromString()` methods
- [ ] Run `DateTime_test.cpp`
- [ ] Commit: "Refactor: Replace QDateTime with std::chrono"

## Phase 2: File System Operations

### File Class
- [ ] Replace `QFile` with `std::fstream`
- [ ] Replace `QFileInfo` with `std::filesystem::path`
- [ ] Replace `QDir` with `std::filesystem`
- [ ] Update `exists()` method
- [ ] Update `readable()` method
- [ ] Update `writable()` method
- [ ] Update `isDirectory()` method
- [ ] Update `fileSize()` method
- [ ] Update `find()` method
- [ ] Update `rename()` method
- [ ] Run `File_test.cpp`
- [ ] Commit: "Refactor: Replace Qt file operations with std::filesystem"

### DocumentIdentifier
- [ ] Update path handling to use std::filesystem
- [ ] Remove QDir usage
- [ ] Run tests
- [ ] Commit: "Refactor: Remove QDir from DocumentIdentifier"

## Phase 3: Process Management

### ExternalProcess
- [ ] Replace QProcess with boost::process or custom wrapper
- [ ] Update `run()` method signature (QString → std::string)
- [ ] Update process output reading
- [ ] Update environment variable handling
- [ ] Run `ExternalProcess_test.cpp`
- [ ] Commit: "Refactor: Replace QProcess with boost::process"

### RWrapper
- [ ] Update to use new ExternalProcess API
- [ ] Replace QString/QStringList parameters
- [ ] Run R integration tests
- [ ] Commit: "Refactor: Remove Qt from RWrapper"

### PythonInfo
- [ ] Replace QProcess with new API
- [ ] Replace QDir::isRelativePath
- [ ] Replace QStringList
- [ ] Run Python integration tests
- [ ] Commit: "Refactor: Remove Qt from PythonInfo"

### JavaInfo
- [ ] Replace QProcess with new API
- [ ] Replace QDir usage
- [ ] Run Java integration tests
- [ ] Commit: "Refactor: Remove Qt from JavaInfo"

## Phase 4: Networking

### NetworkGetRequest
- [ ] Check status of PR #8201
- [ ] If not done: Replace QNetworkAccessManager with cpp-httplib
- [ ] Update `getResponse()` return type
- [ ] Update `getResponseBinary()` return type
- [ ] Run network tests
- [ ] Commit: "Refactor: Replace Qt networking with cpp-httplib"

### UpdateCheck
- [ ] Update to use new NetworkGetRequest API
- [ ] Replace QUrl, QFile, QDir
- [ ] Replace QDateTime for file modification time
- [ ] Run update check tests
- [ ] Commit: "Refactor: Remove Qt from UpdateCheck"

### MascotRemoteQuery
- [ ] Update networking code
- [ ] Replace QRegularExpression with std::regex
- [ ] Run Mascot tests
- [ ] Commit: "Refactor: Remove Qt from MascotRemoteQuery"

## Phase 5: File Watching

### FileWatcher
- [ ] Implement platform-specific file watching
  - [ ] Linux: inotify
  - [ ] Windows: ReadDirectoryChangesW
  - [ ] macOS: FSEvents
- [ ] Or integrate cross-platform library (efsw)
- [ ] Update API to remove QFileSystemWatcher
- [ ] Run FileWatcher tests
- [ ] Commit: "Refactor: Replace QFileSystemWatcher with platform APIs"

## Phase 6: Application Framework

### TOPPBase
- [ ] Replace QString in method signatures
- [ ] Update parameter handling
- [ ] Run TOPP tool tests
- [ ] Commit: "Refactor: Remove Qt from TOPPBase"

### ToolHandler
- [ ] Remove QStringList includes
- [ ] Update tool listing methods
- [ ] Run tests
- [ ] Commit: "Refactor: Remove Qt from ToolHandler"

### ConsoleUtils
- [ ] Update string handling
- [ ] Run tests
- [ ] Commit: "Refactor: Remove Qt from ConsoleUtils"

## Phase 7: Format Handlers & Misc

### ExperimentalDesign
- [ ] Replace QString usage
- [ ] Remove QFileInfo includes
- [ ] Run tests
- [ ] Commit: "Refactor: Remove Qt from ExperimentalDesign"

### FIAMSScheduler
- [ ] Remove QDir includes
- [ ] Update path handling
- [ ] Run tests
- [ ] Commit: "Refactor: Remove Qt from FIAMSScheduler"

### Other Format Files
- [ ] Update remaining FORMAT files (15+ files)
- [ ] Run format tests
- [ ] Commit: "Refactor: Remove Qt from format handlers"

## Phase 8: CMake Configuration

### OpenMS CMakeLists.txt
- [ ] Remove `Qt6::Core` from `OPENMS_DEP_LIBRARIES`
- [ ] Remove `Qt6::Network` from `OPENMS_DEP_LIBRARIES`
- [ ] Verify find_package(Qt6) is still called for GUI
- [ ] Update any Qt-specific compile definitions
- [ ] Commit: "Build: Remove Qt dependencies from OpenMS core library"

### Test CMake
- [ ] Build project with `-DWITH_GUI=OFF`
- [ ] Verify build succeeds without Qt
- [ ] Build project with `-DWITH_GUI=ON`
- [ ] Verify GUI still builds with Qt

## Phase 9: GUI Boundary Layer

### Qt Conversion Utilities
- [ ] Create `QtConversion.h` header (optional)
- [ ] Add conversion helpers for openms_gui
- [ ] Update openms_gui files to use helpers
- [ ] Test GUI functionality
- [ ] Commit: "GUI: Add Qt conversion utilities for boundary layer"

## Phase 10: Testing & Validation

### Build Tests
- [ ] Clean build directory
- [ ] Build with `-DWITH_GUI=OFF`: `cmake --build build --target OpenMS`
- [ ] Build with `-DWITH_GUI=ON`: `cmake --build build`
- [ ] Verify no Qt linking in core library: `ldd build/lib/libOpenMS.so | grep -i qt`
- [ ] Verify Qt linking in GUI library: `ldd build/lib/libOpenMS_GUI.so | grep -i qt`

### Unit Tests
- [ ] Run full test suite: `ctest --test-dir build -j4`
- [ ] Fix any failing tests
- [ ] Verify test coverage hasn't decreased

### Python Bindings
- [ ] Build pyOpenMS
- [ ] Run pyOpenMS tests
- [ ] Verify import works: `python -c "import pyopenms"`

### GUI Tests
- [ ] Launch TOPPView
- [ ] Test basic GUI functionality
- [ ] Verify Qt widgets work correctly

### Performance Testing
- [ ] Measure build time (should be < 2 hours with vcpkg)
- [ ] Run performance benchmarks
- [ ] Compare with baseline

### Documentation
- [ ] Update CHANGELOG
- [ ] Update any Qt-related documentation
- [ ] Update developer documentation if needed

## Phase 11: Pull Request

### Prepare PR
- [ ] Rebase on latest develop branch
- [ ] Squash commits if needed (or keep logical commits)
- [ ] Write comprehensive PR description
- [ ] Add "Fixes #8200" to PR body
- [ ] List all major changes
- [ ] Add testing notes

### PR Checklist Items
- [ ] Code follows OpenMS coding conventions
- [ ] All tests pass
- [ ] No memory leaks (valgrind clean)
- [ ] Documentation updated
- [ ] CHANGELOG updated
- [ ] CMake files updated
- [ ] GUI functionality preserved
- [ ] pyOpenMS builds successfully

### Submit PR
- [ ] Push branch to fork: `git push origin feature/remove-qt-from-core`
- [ ] Create PR on GitHub
- [ ] Tag maintainers: @jpfeuffer @timosachsenberg
- [ ] Respond to review comments
- [ ] Make requested changes
- [ ] Get approval
- [ ] Merge!

## Success Metrics

### Quantitative
- [ ] Build time reduced (target: < 2 hours)
- [ ] 0 Qt includes in src/openms/ (excluding extern)
- [ ] 0 Qt includes in src/openswathalgo/
- [ ] 0 Qt dependencies in OpenMS CMake
- [ ] 100% test pass rate maintained
- [ ] pyOpenMS package size reduced

### Qualitative
- [ ] Code is more maintainable
- [ ] Clearer separation between core and GUI
- [ ] Easier to deploy Python packages
- [ ] No Qt license concerns for core library
- [ ] Easier vcpkg builds on Windows

---

## Notes

### Important Reminders
- ⚠️ **openms_gui must keep Qt dependencies**
- ⚠️ Test after each phase, don't accumulate changes
- ⚠️ Keep commits small and focused
- ⚠️ Check PR #8201 status before working on networking

### Helpful Commands
```bash
# Build core without GUI
cmake -B build -S . -DWITH_GUI=OFF -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target OpenMS -j4

# Run specific test
ctest --test-dir build -R String_test -V

# Check Qt linkage
ldd build/lib/libOpenMS.so | grep -i qt

# Run audit script
bash tools/scripts/audit_qt_dependencies.sh
```

### Resources
- [Issue #8200](https://github.com/OpenMS/OpenMS/issues/8200)
- [PR #8201 - Networking](https://github.com/OpenMS/OpenMS/pull/8201)
- [OpenMS Developer Docs](https://openms.de/documentation/)
- [C++ Filesystem Reference](https://en.cppreference.com/w/cpp/filesystem)
- [Boost.Process](https://www.boost.org/doc/libs/1_83_0/doc/html/process.html)
