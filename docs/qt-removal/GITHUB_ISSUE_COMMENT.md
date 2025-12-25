# GitHub Issue Comment - Ready to Work on Issue #8200

Hi @jpfeuffer and @timosachsenberg,

I'd like to work on this issue! I've set up my development environment and conducted a thorough audit of the codebase to understand the scope.

## Audit Summary

I've analyzed the Qt dependencies across the OpenMS core library:

**Statistics:**
- **60 files** require Qt removal
- **Major Qt usage:**
  - QString: 134 occurrences
  - QDir: 61 occurrences  
  - QStringList: 48 occurrences
  - QDateTime: 30 occurrences
  - QProcess: 29 occurrences
  - QFileInfo: 26 occurrences
  - QByteArray: 25 occurrences
  - QFile: 18 occurrences
  - QNetworkAccessManager: 7 occurrences

**Good news:** OpenSwathAlgo is already Qt-free! ✅

## Proposed Implementation Plan

I've broken this down into manageable phases:

### Phase 1: Core Data Structures (Week 1-2)
- Refactor `String`, `DataValue`, `StringUtils` classes
- Replace `QDateTime` with `std::chrono`
- Remove `toQString()` methods and QString constructors

### Phase 2: File System Operations (Week 2-3)
- Replace `QFile`, `QFileInfo`, `QDir` with `std::filesystem` (C++17)
- Update `File.cpp` and related classes
- Implement equivalent functionality using standard library

### Phase 3: Process Management (Week 3-4)
- Replace `QProcess` in `ExternalProcess` with `boost::process` (already used in OpenMS)
- Update `RWrapper`, `PythonInfo`, `JavaInfo` to use new API

### Phase 4: Networking (Week 4-5)
- Coordinate with PR #8201 which addresses `QNetworkAccessManager`
- Use `cpp-httplib` or existing solution from that PR
- Update `NetworkGetRequest` and `UpdateCheck`

### Phase 5: File Watching (Week 5)
- Replace `QFileSystemWatcher` with platform-specific APIs (inotify on Linux)
- Consider cross-platform library like `efsw` or `boost::asio`

### Phase 6-7: Application Framework & CMake (Week 6)
- Update `TOPPBase`, `ToolHandler`
- Remove `Qt6::Core` and `Qt6::Network` from `src/openms/CMakeLists.txt`

### Phase 8: GUI Boundary Layer (Week 7)
- Create conversion utilities for `openms_gui` boundary
- Ensure GUI library still functions with Qt

## Development Approach

- **Small, focused commits**: Each phase will be a separate commit or small PR
- **Comprehensive testing**: Run unit tests after each change
- **No functionality regression**: All existing behavior preserved
- **Documentation**: Update as needed

## Expected Benefits

After completion:
- ✅ Faster build times (< 2 hours with vcpkg)
- ✅ Smaller pyOpenMS packages (no Qt dependencies)
- ✅ No Qt licensing concerns for core library
- ✅ Cleaner separation between core and GUI
- ✅ Easier deployment for headless servers

## Questions for Maintainers

1. **boost::process**: Is this acceptable for replacing QProcess, or would you prefer a different approach?
2. **std::filesystem**: C++17 filesystem is used - is this OK, or do we need C++14 compatibility?
3. **PR #8201**: Should I wait for the networking PR to merge, or coordinate/integrate that work?
4. **PR Strategy**: Would you prefer one large PR or multiple smaller PRs for each phase?
5. **Testing priority**: Are there specific tools/workflows that are critical to test thoroughly?

## My Setup

- OS: Linux (Ubuntu dual-boot)
- Development: VS Code with C++/CMake tools
- Experience: Proficient in C++, familiar with standard library and Boost
- Already cloned and built OpenMS successfully

I'm ready to start as soon as I get feedback on the approach. Looking forward to contributing to OpenMS!

Best regards,
Prachi
