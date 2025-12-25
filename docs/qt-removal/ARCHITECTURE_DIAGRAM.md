# OpenMS Architecture & Qt Dependency Map

## Current Architecture (Before)

```
┌─────────────────────────────────────────────────────────────┐
│                     OpenMS Project                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌────────────────────────────────────────────────────┐   │
│  │  src/openms/  (Core Library)                       │   │
│  │  ❌ Qt6::Core + Qt6::Network                       │   │
│  │                                                     │   │
│  │  Contains:                                          │   │
│  │  • DATASTRUCTURES/String → Uses QString            │   │
│  │  • DATASTRUCTURES/DataValue → Uses QString         │   │
│  │  • SYSTEM/File → Uses QFile, QDir, QFileInfo      │   │
│  │  • SYSTEM/ExternalProcess → Uses QProcess          │   │
│  │  • SYSTEM/NetworkGetRequest → Uses QNetwork        │   │
│  │  • FORMAT/* → Various Qt usage                     │   │
│  │                                                     │   │
│  │  Dependencies: Boost, XercesC, LibSVM, Qt6 ❌      │   │
│  └────────────────────────────────────────────────────┘   │
│                          ↑                                  │
│                          │ Links to                         │
│                          │                                  │
│  ┌────────────────────────────────────────────────────┐   │
│  │  src/openswathalgo/  (OpenSwath Library)           │   │
│  │  ✅ Already Qt-free                                │   │
│  └────────────────────────────────────────────────────┘   │
│                          ↑                                  │
│                          │ Used by                          │
│  ┌───────────────────┬──┴─────────────┬─────────────┐     │
│  │                   │                │             │     │
│  │  src/topp/        │  src/pyOpenMS/ │ src/openms_gui/ │
│  │  (CLI Tools)      │  (Python)      │ (GUI Library)   │
│  │  ❌ Needs Qt      │  ❌ Pulls Qt   │ ✅ Keep Qt      │
│  │  currently        │  unnecessarily │ (GUI needs it)  │
│  └───────────────────┴────────────────┴─────────────┘     │
│                                                             │
└─────────────────────────────────────────────────────────────┘

Problems:
• pyOpenMS users forced to install Qt
• TOPP CLI tools require Qt for no GUI reason
• Slower build times
• Qt licensing concerns
• Harder to deploy headless
```

## Target Architecture (After Issue #8200)

```
┌─────────────────────────────────────────────────────────────┐
│                     OpenMS Project                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌────────────────────────────────────────────────────┐   │
│  │  src/openms/  (Core Library)                       │   │
│  │  ✅ NO Qt Dependencies                             │   │
│  │                                                     │   │
│  │  Contains:                                          │   │
│  │  • DATASTRUCTURES/String → std::string             │   │
│  │  • DATASTRUCTURES/DataValue → std::string          │   │
│  │  • SYSTEM/File → std::filesystem                   │   │
│  │  • SYSTEM/ExternalProcess → boost::process         │   │
│  │  • SYSTEM/NetworkGetRequest → cpp-httplib          │   │
│  │  • FORMAT/* → Standard C++ only                    │   │
│  │                                                     │   │
│  │  Dependencies: Boost, XercesC, LibSVM ✅           │   │
│  └────────────────────────────────────────────────────┘   │
│                          ↑                                  │
│                          │ Links to                         │
│                          │                                  │
│  ┌────────────────────────────────────────────────────┐   │
│  │  src/openswathalgo/  (OpenSwath Library)           │   │
│  │  ✅ Already Qt-free                                │   │
│  └────────────────────────────────────────────────────┘   │
│                          ↑                                  │
│                          │ Used by                          │
│  ┌───────────────────┬──┴─────────────┬─────────────┐     │
│  │                   │                │             │     │
│  │  src/topp/        │  src/pyOpenMS/ │ src/openms_gui/ │
│  │  (CLI Tools)      │  (Python)      │ (GUI Library)   │
│  │  ✅ No Qt needed  │  ✅ No Qt      │ ✅ Keep Qt      │
│  │  Clean CLI        │  dependency    │ (Uses Qt6)      │
│  │                   │                │                 │
│  │                   │                │ • QString/Qt    │
│  │                   │                │   conversion at │
│  │                   │                │   boundary only │
│  └───────────────────┴────────────────┴─────────────┘     │
│                                                             │
└─────────────────────────────────────────────────────────────┘

Benefits:
✅ pyOpenMS packages 50-80% smaller
✅ TOPP tools are pure CLI (no Qt needed)
✅ Faster builds (< 2 hours with vcpkg)
✅ No Qt licensing issues for core
✅ Easier headless deployment
✅ GUI still fully functional
```

## Key Files & Qt Usage Map

### 🔴 High Priority (Core Data Types - Wide Impact)

```
src/openms/include/OpenMS/DATASTRUCTURES/
├── String.h ⚠️ QString                    [PHASE 1 - Week 1]
├── DataValue.h ⚠️ QString                 [PHASE 1 - Week 1]
├── StringUtils.h ⚠️ QString conversions   [PHASE 1 - Week 1]
└── StringListUtils.h ⚠️ QStringList       [PHASE 1 - Week 1]

src/openms/source/DATASTRUCTURES/
├── String.cpp ⚠️ toQString(), QString()
├── DataValue.cpp ⚠️ toQString(), QString()
├── DateTime.cpp ⚠️ QDateTime               [PHASE 1 - Week 2]
└── StringListUtils.cpp ⚠️ fromQStringList()
```

### 🟡 Medium Priority (System Utilities)

```
src/openms/include/OpenMS/SYSTEM/
├── File.h ⚠️ QFile, QDir, QFileInfo       [PHASE 2 - Week 2-3]
├── ExternalProcess.h ⚠️ QProcess          [PHASE 3 - Week 3-4]
├── NetworkGetRequest.h ⚠️ QNetwork        [PHASE 4 - Week 4-5]
├── FileWatcher.h ⚠️ QFileSystemWatcher    [PHASE 5 - Week 5]
├── RWrapper.h ⚠️ QProcess, QString
├── UpdateCheck.cpp ⚠️ Multiple Qt types
├── PythonInfo.cpp ⚠️ QProcess, QDir
└── JavaInfo.cpp ⚠️ QProcess, QDir

src/openms/source/SYSTEM/
├── File.cpp ⚠️ 18 QFile usages
├── ExternalProcess.cpp ⚠️ 29 QProcess usages
├── NetworkGetRequest.cpp ⚠️ 7 QNetworkAccessManager
└── FileWatcher.cpp ⚠️ 7 QFileSystemWatcher
```

### 🟢 Lower Priority (Application Layer)

```
src/openms/include/OpenMS/APPLICATIONS/
├── TOPPBase.h ⚠️ QString                  [PHASE 6 - Week 6]
└── ToolHandler.h ⚠️ QStringList

src/openms/source/APPLICATIONS/
├── TOPPBase.cpp
├── ToolHandler.cpp
└── ConsoleUtils.cpp

src/openms/source/FORMAT/
├── MascotRemoteQuery.cpp ⚠️ QRegularExpression
├── Various other format files (~15 files)
└── ...
```

### ⚪ Untouched (Keep Qt)

```
src/openms_gui/                            [NO CHANGES]
├── ✅ Keep all Qt dependencies
├── ✅ Still uses QString, QWidget, etc.
└── ✅ Add conversion helpers at boundary
```

## Replacement Strategy Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                  Qt Class Replacements                      │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  QString                                                    │
│  ├─ std::string                        (C++11)             │
│  └─ For UTF-16: std::u16string         (rare case)         │
│                                                             │
│  QStringList                                                │
│  └─ std::vector<std::string>           (C++11)             │
│                                                             │
│  QFile                                                      │
│  ├─ std::ifstream / std::ofstream      (C++98)             │
│  └─ std::filesystem::path operations   (C++17)             │
│                                                             │
│  QFileInfo                                                  │
│  └─ std::filesystem::path + operations (C++17)             │
│                                                             │
│  QDir                                                       │
│  ├─ std::filesystem::directory_iterator (C++17)            │
│  ├─ std::filesystem::create_directories (C++17)            │
│  └─ std::filesystem::path operations    (C++17)            │
│                                                             │
│  QDateTime                                                  │
│  ├─ std::chrono::system_clock          (C++11)             │
│  ├─ std::format (for formatting)       (C++20)             │
│  └─ OR boost::date_time                (if complex)        │
│                                                             │
│  QByteArray                                                 │
│  ├─ std::vector<char>                  (binary data)       │
│  └─ std::string                        (text data)         │
│                                                             │
│  QProcess                                                   │
│  ├─ boost::process                     (recommended)       │
│  ├─ fork() + exec()                    (Unix only)         │
│  └─ CreateProcess()                    (Windows only)      │
│                                                             │
│  QList<T> / QVector<T>                                      │
│  └─ std::vector<T>                     (C++98)             │
│                                                             │
│  QMap<K, V>                                                 │
│  ├─ std::map<K, V>                     (ordered)           │
│  └─ std::unordered_map<K, V>           (hash-based)        │
│                                                             │
│  QNetworkAccessManager                                      │
│  ├─ cpp-httplib                        (PR #8201)          │
│  └─ OR boost::asio                     (alternative)       │
│                                                             │
│  QFileSystemWatcher                                         │
│  ├─ inotify API                        (Linux)             │
│  ├─ ReadDirectoryChangesW              (Windows)           │
│  ├─ FSEvents                           (macOS)             │
│  └─ OR efsw library                    (cross-platform)    │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## Implementation Phases Timeline

```
Week 1-2: PHASE 1 - Core Data Structures
│
├─ Day 1-2:   String & StringUtils refactoring
├─ Day 3-4:   DataValue refactoring  
├─ Day 5-7:   DateTime refactoring
└─ Day 8-10:  Testing & bug fixes
              └─ Commit: "Refactor: Remove Qt from core data types"

Week 2-3: PHASE 2 - File System
│
├─ Day 1-3:   File.cpp QFile → std::fstream
├─ Day 4-6:   File.cpp QDir/QFileInfo → std::filesystem
├─ Day 7-8:   DocumentIdentifier & other file users
└─ Day 9-10:  Testing
              └─ Commit: "Refactor: Replace Qt file I/O with std::filesystem"

Week 3-4: PHASE 3 - Process Management
│
├─ Day 1-4:   ExternalProcess QProcess → boost::process
├─ Day 5-6:   RWrapper updates
├─ Day 7-8:   PythonInfo, JavaInfo updates
└─ Day 9-10:  Testing
              └─ Commit: "Refactor: Replace QProcess with boost::process"

Week 4-5: PHASE 4 - Networking
│
├─ Day 1-2:   Check PR #8201 status
├─ Day 3-6:   NetworkGetRequest refactoring (if needed)
├─ Day 7-8:   UpdateCheck, MascotRemoteQuery updates
└─ Day 9-10:  Testing
              └─ Commit: "Refactor: Remove Qt networking dependencies"

Week 5: PHASE 5 - File Watching
│
├─ Day 1-4:   FileWatcher platform implementations
├─ Day 5-7:   Integration & testing
└─ Day 8-10:  Edge case handling
              └─ Commit: "Refactor: Replace QFileSystemWatcher"

Week 6: PHASE 6-7 - Application & CMake
│
├─ Day 1-3:   TOPPBase refactoring
├─ Day 4-5:   ToolHandler, ConsoleUtils
├─ Day 6-7:   Remaining FORMAT files
├─ Day 8-9:   CMakeLists.txt updates
└─ Day 10:    Build testing
              └─ Commit: "Build: Remove Qt from OpenMS CMakeLists"

Week 7-8: PHASE 8-9 - GUI Boundary & Final Testing
│
├─ Day 1-3:   Create QtConversion.h helpers
├─ Day 4-5:   Update openms_gui boundary code
├─ Day 6-7:   Full test suite
├─ Day 8-9:   Performance testing & optimization
├─ Day 10-12: pyOpenMS testing
└─ Day 13-14: Documentation & PR preparation
              └─ PR: "Remove all Qt dependencies from OpenMS core"
```

## Testing Strategy

```
┌─────────────────────────────────────────────────────────────┐
│                     Testing Pyramid                         │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│                      ▲                                      │
│                     ╱ ╲  Integration Tests                  │
│                    ╱   ╲  • Full build                      │
│                   ╱     ╲ • pyOpenMS                        │
│                  ╱       ╲• GUI functionality               │
│                 ╱─────────╲                                 │
│                ╱           ╲                                │
│               ╱  Unit Tests ╲                               │
│              ╱  • String_test ╲                             │
│             ╱   • File_test    ╲                            │
│            ╱    • Process_test  ╲                           │
│           ╱      • After each    ╲                          │
│          ╱        major change    ╲                         │
│         ╱───────────────────────────╲                       │
│        ╱                             ╲                      │
│       ╱     Compile-Time Checks       ╲                     │
│      ╱    • Build after each commit    ╲                    │
│     ╱     • No Qt symbols in libOpenMS  ╲                   │
│    ╱      • Linker verification          ╲                  │
│   ╱────────────────────────────────────────╲                │
│                                                             │
└─────────────────────────────────────────────────────────────┘

Test Commands:
1. After each file change:    cmake --build build --target OpenMS
2. After logical unit:        ctest -R SpecificTest_test
3. After each phase:          ctest --test-dir build
4. Before commit:             Full build + full tests
5. Before PR:                 Build with GUI + pyOpenMS
```

## Success Metrics

```
┌─────────────────────────────────────────────────────────────┐
│                   Before → After                            │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Build Time (with vcpkg on Windows)                        │
│  4-6 hours  ───────────────────→  < 2 hours    ✅          │
│                                                             │
│  pyOpenMS Package Size                                      │
│  ~300 MB    ───────────────────→  ~100 MB      ✅          │
│                                                             │
│  Core Library Dependencies                                  │
│  Qt6::Core, Qt6::Network  ─────→  None          ✅          │
│                                                             │
│  Files with Qt in src/openms/                               │
│  60 files   ───────────────────→  0 files       ✅          │
│                                                             │
│  QString Occurrences                                        │
│  134        ───────────────────→  0             ✅          │
│                                                             │
│  Qt License Concerns                                        │
│  LGPL       ───────────────────→  No Qt         ✅          │
│                                                             │
│  Test Pass Rate                                             │
│  100%       ───────────────────→  100%          ✅          │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

---

## Quick Reference: What Goes Where

### ❌ Remove Qt From (Phase 1-7)
- `src/openms/include/OpenMS/**/*.h`
- `src/openms/source/**/*.cpp`
- `src/openms/CMakeLists.txt`

### ✅ Keep Qt In (Untouched)
- `src/openms_gui/**/*` (All GUI code)
- TOPP GUI tools (TOPPView, etc.)

### 🔄 Optional Boundary Helpers
- New: `src/openms/include/OpenMS/DATASTRUCTURES/QtConversion.h`
  - Only compiled with `-DWITH_GUI=ON`
  - Provides std::string ↔ QString helpers for GUI

---

**This diagram provides a visual overview of the entire refactoring task. Refer to it when planning your work!**
