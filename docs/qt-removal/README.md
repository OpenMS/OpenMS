# VS Code Setup for OpenMS Qt Removal (Issue #8200)

This directory contains all the resources you need to work on removing Qt dependencies from the OpenMS core library.

## 📁 Files in This Directory

### Configuration Files
- **`settings.json`** - VS Code C++/CMake settings optimized for OpenMS development
- **`tasks.json`** - Build tasks for compiling OpenMS without GUI and running tests

### Documentation
- **`QT_REMOVAL_ACTION_PLAN.md`** - Comprehensive 8-phase implementation plan
- **`QT_REMOVAL_CHECKLIST.md`** - Detailed checklist with checkboxes to track progress
- **`QT_REPLACEMENT_PATTERNS.md`** - Code patterns for Qt → Standard C++ replacements
- **`GITHUB_ISSUE_COMMENT.md`** - Draft comment to post on GitHub Issue #8200

### Scripts
- **`../tools/scripts/audit_qt_dependencies.sh`** - Audit script to find all Qt usage

## 🚀 Quick Start

### 1. Initial Setup

```bash
# Navigate to workspace
cd /home/prachi/openMS/OpenMS

# Create your feature branch
git checkout develop
git pull upstream develop
git checkout -b feature/remove-qt-from-core

# Run the Qt audit
bash tools/scripts/audit_qt_dependencies.sh
```

### 2. Build Baseline

First, make sure you can build the project successfully:

```bash
# Configure (without GUI for faster builds during development)
cmake -B build -S . -DWITH_GUI=OFF -DCMAKE_BUILD_TYPE=Debug

# Build OpenMS library
cmake --build build --target OpenMS -j4

# Run tests to establish baseline
ctest --test-dir build -j4
```

### 3. Start Development

1. Open VS Code: `code .`
2. Review the action plan: `.vscode/QT_REMOVAL_ACTION_PLAN.md`
3. Use the checklist: `.vscode/QT_REMOVAL_CHECKLIST.md`
4. Reference replacement patterns: `.vscode/QT_REPLACEMENT_PATTERNS.md`

### 4. Post on GitHub

Copy the content from `.vscode/GITHUB_ISSUE_COMMENT.md` and post it on [Issue #8200](https://github.com/OpenMS/OpenMS/issues/8200) to let maintainers know you're working on it.

## 🔧 VS Code Tasks

Press `Ctrl+Shift+B` (or `Cmd+Shift+B` on Mac) to see available build tasks:

- **CMake: Configure (No GUI)** - Configure build without GUI components
- **CMake: Build OpenMS Library** - Build just the OpenMS core library
- **Run Qt Dependency Audit** - Re-run the Qt audit script
- **Clean Build Directory** - Clean the build directory

## 📊 Current Status (From Audit)

**Total Files Needing Changes:** 60

**Qt Class Usage:**
- QString: 134 occurrences
- QDir: 61 occurrences
- QStringList: 48 occurrences
- QDateTime: 30 occurrences
- QProcess: 29 occurrences
- QFileInfo: 26 occurrences
- QByteArray: 25 occurrences
- QFile: 18 occurrences
- QNetworkAccessManager: 7 occurrences
- QFileSystemWatcher: 7 occurrences

**✅ Good News:** OpenSwathAlgo is already Qt-free!

## 📖 Implementation Strategy

### Recommended Order

1. **Start Small**: Begin with `String` and `DataValue` classes (high impact, well-contained)
2. **Test Frequently**: Run unit tests after each major change
3. **Commit Often**: Make small, focused commits for each logical change
4. **Phase-by-Phase**: Follow the 8 phases in the action plan
5. **Document**: Update CHANGELOG and comments as you go

### Example Workflow

```bash
# 1. Make changes to String class
vim src/openms/include/OpenMS/DATASTRUCTURES/String.h
vim src/openms/source/DATASTRUCTURES/String.cpp

# 2. Build
cmake --build build --target OpenMS -j4

# 3. Test
ctest --test-dir build -R String_test -V

# 4. Commit if tests pass
git add src/openms/include/OpenMS/DATASTRUCTURES/String.h
git add src/openms/source/DATASTRUCTURES/String.cpp
git commit -m "Refactor: Remove QString from String class"

# 5. Re-run audit to track progress
bash tools/scripts/audit_qt_dependencies.sh
```

## 🔍 Helpful VS Code Shortcuts

- `Ctrl+Shift+F` - Global search (find all Qt usage)
- `Ctrl+H` - Find and replace in current file
- `Ctrl+Shift+H` - Global find and replace
- `F12` - Go to definition
- `Shift+F12` - Find all references
- `Ctrl+P` - Quick file open
- `Ctrl+`` - Toggle terminal

## 🧪 Testing Strategy

### After Each Change

1. **Build**: `cmake --build build --target OpenMS -j4`
2. **Test**: `ctest --test-dir build -R <TestName> -V`
3. **Check Errors**: Use VS Code's Problems panel or `get_errors` tool

### Before Committing

1. **Full Build**: Build entire project including tests
2. **Full Test Suite**: Run all tests
3. **Check Qt Linkage**: `ldd build/lib/libOpenMS.so | grep -i qt` (should be empty)
4. **Memory Check** (optional): `valgrind --leak-check=full ./test_executable`

### Before Creating PR

1. **Build with GUI**: `cmake -B build -S . -DWITH_GUI=ON && cmake --build build`
2. **GUI Tests**: Verify openms_gui still works
3. **pyOpenMS**: Test Python bindings if possible
4. **Documentation**: Update CHANGELOG

## 📚 Key Resources

### Documentation
- [Issue #8200](https://github.com/OpenMS/OpenMS/issues/8200)
- [Related PR #8201 - Networking](https://github.com/OpenMS/OpenMS/pull/8201)
- [OpenMS Developer Docs](https://openms.de/documentation/)
- [OpenMS Contributing Guide](../CONTRIBUTING.md)

### C++ References
- [std::filesystem Reference](https://en.cppreference.com/w/cpp/filesystem)
- [std::chrono Reference](https://en.cppreference.com/w/cpp/chrono)
- [Boost.Process Docs](https://www.boost.org/doc/libs/1_83_0/doc/html/process.html)
- [Boost.DateTime Docs](https://www.boost.org/doc/libs/1_83_0/doc/html/date_time.html)

## ⚠️ Important Notes

### DO
- ✅ Remove Qt from `src/openms/` (core library)
- ✅ Remove Qt from `src/openswathalgo/` (already done!)
- ✅ Test after every change
- ✅ Keep commits small and focused
- ✅ Document your changes

### DON'T
- ❌ Touch `src/openms_gui/` - GUI must keep Qt!
- ❌ Break existing functionality
- ❌ Make huge commits with 50+ file changes
- ❌ Skip testing
- ❌ Forget to update CMakeLists.txt

## 🎯 Success Criteria

Your PR is ready when:

- [ ] No Qt includes in `src/openms/` (except extern)
- [ ] No Qt includes in `src/openswathalgo/`
- [ ] `src/openms/CMakeLists.txt` doesn't link Qt
- [ ] All unit tests pass
- [ ] pyOpenMS builds successfully
- [ ] GUI library and tools still work
- [ ] No functionality regression
- [ ] Build time improved
- [ ] CHANGELOG updated

## 🆘 Getting Help

### If You're Stuck
1. Review `.vscode/QT_REPLACEMENT_PATTERNS.md` for common patterns
2. Search existing code for similar patterns
3. Check Boost documentation
4. Ask on GitHub issue or OpenMS Gitter/Discord

### Common Issues
- **Compile errors after removing Qt**: Check for implicit QString conversions
- **Tests failing**: Ensure string formatting/parsing matches original behavior
- **Linker errors**: Update CMakeLists.txt dependencies
- **Performance regression**: Profile and optimize critical paths

## 📞 Contact

- **GitHub**: Comment on [Issue #8200](https://github.com/OpenMS/OpenMS/issues/8200)
- **Maintainers**: @jpfeuffer, @timosachsenberg
- **Your GitHub**: @Prachi194agrawal

---

## Quick Commands Cheat Sheet

```bash
# Development
code .                                    # Open in VS Code
bash tools/scripts/audit_qt_dependencies.sh  # Run audit

# Building
cmake -B build -S . -DWITH_GUI=OFF       # Configure
cmake --build build --target OpenMS -j4  # Build
ctest --test-dir build -R Test_test -V   # Run test

# Git
git status                               # Check status
git add <files>                          # Stage changes
git commit -m "message"                  # Commit
git push origin feature/remove-qt-from-core  # Push to fork

# Verification
ldd build/lib/libOpenMS.so | grep -i qt  # Check Qt linkage (should be empty)
nm -D build/lib/libOpenMS.so | grep -i qt  # Check Qt symbols

# Cleanup
rm -rf build                             # Clean build
git clean -fdx                           # Clean all untracked files (careful!)
```

---

**Good luck with your contribution to OpenMS! 🚀**
