# Step-by-Step Guide: Creating Your Pull Request

This guide will walk you through creating and submitting your pull request for Issue #8200.

## Current Status ✅

- ✅ Feature branch created: `feature/remove-qt-from-core`
- ✅ Fork configured as `origin`: https://github.com/Prachi194agrawal/OpenMS
- ✅ Upstream configured: https://github.com/OpenMS/OpenMS.git
- ✅ VS Code workspace set up
- ✅ Documentation and tools created

## Step 1: Commit Your Setup Files

First, let's commit all the setup work you've done:

```bash
cd /home/prachi/openMS/OpenMS

# Add all the setup files
git add .vscode/
git add tools/scripts/audit_qt_dependencies.sh

# Check what will be committed
git status

# Commit with a descriptive message
git commit -m "Setup: Add VS Code configuration and documentation for Qt removal (Issue #8200)

- Add comprehensive action plan with 8 implementation phases
- Add detailed checklist for tracking progress
- Add Qt replacement patterns and code examples
- Add architecture diagrams and visual guides
- Add audit script for tracking Qt dependencies
- Add VS Code tasks for building and testing
- Configure C++/CMake settings for development

This setup prepares the development environment for removing Qt
dependencies from OpenMS core library as described in Issue #8200."
```

## Step 2: Push to Your Fork

```bash
# Push your feature branch to your fork
git push origin feature/remove-qt-from-core
```

## Step 3: Create Pull Request on GitHub

### Option A: Using GitHub CLI (if installed)
```bash
gh pr create \
  --title "WIP: Setup for removing Qt dependencies from OpenMS core (Issue #8200)" \
  --body-file .vscode/GITHUB_ISSUE_COMMENT.md \
  --base develop \
  --head Prachi194agrawal:feature/remove-qt-from-core
```

### Option B: Using GitHub Web Interface (Recommended)

1. **Go to your fork on GitHub:**
   https://github.com/Prachi194agrawal/OpenMS

2. **You'll see a notification:**
   "feature/remove-qt-from-core had recent pushes"
   Click "Compare & pull request"

3. **Fill in the PR details:**

   **Title:**
   ```
   WIP: Setup for removing Qt dependencies from OpenMS core (Issue #8200)
   ```

   **Description:** (Copy from `.vscode/GITHUB_ISSUE_COMMENT.md` and modify)
   ```markdown
   ## Description
   
   This PR is the first step in addressing Issue #8200: Remove all Qt dependencies from the OpenMS and OpenSwath libraries.
   
   ### What This PR Contains
   
   This initial PR sets up the development environment and documentation:
   
   - 📋 Comprehensive 8-phase action plan
   - ✅ Detailed checklist for tracking progress
   - 📖 Qt → C++ replacement patterns and examples
   - 🏗️ Architecture diagrams showing current vs. target state
   - 🔍 Audit script to track Qt dependencies
   - ⚙️ VS Code configuration for C++/CMake development
   - 🛠️ Build tasks for development workflow
   
   ### Audit Results
   
   Current Qt usage in OpenMS core:
   - **60 files** need Qt removal
   - **134 QString** occurrences
   - **61 QDir** occurrences
   - **48 QStringList** occurrences
   - **30 QDateTime** occurrences
   - **29 QProcess** occurrences
   
   ✅ **Good news:** OpenSwathAlgo is already Qt-free!
   
   ## Implementation Plan
   
   I propose to implement this in phases:
   
   ### Phase 1: Core Data Structures (Week 1-2)
   - Refactor `String`, `DataValue`, `StringUtils` classes
   - Replace `QDateTime` with `std::chrono`
   - Remove `toQString()` methods and QString constructors
   
   ### Phase 2: File System Operations (Week 2-3)
   - Replace `QFile`, `QFileInfo`, `QDir` with `std::filesystem` (C++17)
   - Update `File.cpp` and related classes
   
   ### Phase 3: Process Management (Week 3-4)
   - Replace `QProcess` with `boost::process` (already used in OpenMS)
   - Update `ExternalProcess`, `RWrapper`, `PythonInfo`, `JavaInfo`
   
   ### Phase 4: Networking (Week 4-5)
   - Coordinate with PR #8201 which addresses `QNetworkAccessManager`
   - Use `cpp-httplib` or existing solution from that PR
   - Update `NetworkGetRequest` and `UpdateCheck`
   
   ### Phase 5: File Watching (Week 5)
   - Replace `QFileSystemWatcher` with platform-specific APIs or cross-platform library
   
   ### Phase 6-7: Application Framework & CMake (Week 6)
   - Update `TOPPBase`, `ToolHandler`, format handlers
   - Remove `Qt6::Core` and `Qt6::Network` from `src/openms/CMakeLists.txt`
   
   ### Phase 8: GUI Boundary & Testing (Week 7-8)
   - Create conversion utilities for `openms_gui` boundary
   - Full test suite and performance testing
   - Verify pyOpenMS builds successfully
   
   ## Proposed Replacements
   
   - `QString` → `std::string`
   - `QStringList` → `std::vector<std::string>`
   - `QFile/QFileInfo` → `std::filesystem::path` (C++17)
   - `QDir` → `std::filesystem::directory_iterator`
   - `QDateTime` → `std::chrono` or `boost::date_time`
   - `QByteArray` → `std::vector<char>` or `std::string`
   - `QProcess` → `boost::process`
   - `QNetworkAccessManager` → `cpp-httplib` (see PR #8201)
   
   ## Questions for Maintainers
   
   1. **boost::process**: Is this acceptable for replacing QProcess?
   2. **std::filesystem**: C++17 is required - is this OK for OpenMS?
   3. **PR #8201**: Should I wait for the networking PR to merge first?
   4. **PR Strategy**: One large PR or multiple smaller PRs for each phase?
   5. **Testing priority**: Which tools/workflows are most critical to test?
   
   ## Expected Benefits
   
   After completion:
   - ✅ Faster build times (< 2 hours with vcpkg)
   - ✅ Smaller pyOpenMS packages (no Qt dependencies)
   - ✅ No Qt licensing concerns for core library
   - ✅ Cleaner separation between core and GUI
   - ✅ Easier deployment for headless servers
   
   ## Important Notes
   
   - ⚠️ **openms_gui must keep Qt dependencies** - This PR only affects core library
   - ⚠️ All functionality will be preserved - no breaking changes to API behavior
   - ⚠️ Comprehensive testing after each phase
   
   ## Checklist
   
   - [x] Development environment set up
   - [x] Qt dependency audit completed
   - [x] Implementation plan documented
   - [x] VS Code tasks configured
   - [ ] Post plan on GitHub Issue #8200 for feedback
   - [ ] Begin Phase 1 implementation after maintainer approval
   
   Fixes #8200
   
   ---
   
   **Note:** This is marked as WIP (Work In Progress) as this PR currently contains only setup and documentation. Actual code changes will follow after maintainer feedback.
   ```

4. **Base branch:** Make sure it's set to `develop` (not `master`)

5. **Reviewers:** Request review from:
   - @jpfeuffer
   - @timosachsenberg

6. **Labels:** Add if available:
   - `enhancement`
   - `refactoring`
   - `help wanted`

7. **Link to Issue:** In the sidebar, link to Issue #8200

8. **Click:** "Create pull request"

## Step 4: Post Comment on Issue #8200

After creating the PR, go to the issue and add a comment:

https://github.com/OpenMS/OpenMS/issues/8200

**Comment:**
```markdown
Hi @jpfeuffer and @timosachsenberg,

I've started working on this issue! I've created PR #XXXX with setup and documentation.

I've done a thorough audit of the codebase and created a detailed implementation plan. The setup includes:

- 📋 8-phase action plan
- ✅ Detailed checklist (60+ files to modify)
- 📖 Qt → C++ replacement patterns
- 🔍 Automated audit script
- 🏗️ Architecture diagrams

**Key findings:**
- 60 files need Qt removal
- Main usage: QString (134), QDir (61), QStringList (48), QDateTime (30), QProcess (29)
- OpenSwathAlgo is already Qt-free ✅

I'm ready to start implementation, but wanted to get your feedback on the approach first. Please review the PR and let me know:

1. Is the phased approach acceptable?
2. Should I submit one large PR or multiple smaller PRs?
3. Any concerns about using `boost::process` or `std::filesystem`?
4. Should I coordinate with PR #8201 for networking?

Looking forward to your feedback!
```

## Step 5: Wait for Feedback

Before starting actual code changes:
- ⏳ Wait for maintainer feedback on your plan
- 💬 Respond to any questions or concerns
- ✅ Get approval to proceed

## Step 6: Start Implementation (After Approval)

Once you get approval:

```bash
# Make sure you're on the feature branch
git checkout feature/remove-qt-from-core

# Start with Phase 1 - String class
vim src/openms/include/OpenMS/DATASTRUCTURES/String.h

# After making changes:
git add src/openms/include/OpenMS/DATASTRUCTURES/String.h
git add src/openms/source/DATASTRUCTURES/String.cpp
git commit -m "Refactor: Remove QString from String class

- Remove String::toQString() method
- Remove String(const QString&) constructor
- Update all QString-related methods to use std::string
- Add tests to verify functionality preserved

Part of Issue #8200 - Phase 1: Core Data Structures"

# Push updates
git push origin feature/remove-qt-from-core
```

## Alternative: Create Issue Comment First

If you want feedback BEFORE creating a PR, you can:

1. **Just commit locally** (don't push yet)
2. **Post on Issue #8200** with your plan (use `.vscode/GITHUB_ISSUE_COMMENT.md`)
3. **Wait for feedback**
4. **Then push and create PR**

This approach is more conservative and gets buy-in first.

## Troubleshooting

### If push is rejected
```bash
# Pull latest changes from upstream
git fetch upstream
git rebase upstream/develop

# Force push if needed (only on your branch!)
git push origin feature/remove-qt-from-core --force-with-lease
```

### If you need to update your PR description
- Go to the PR on GitHub
- Click "Edit" next to the title
- Update the description
- Click "Save"

## Next Steps After PR Created

1. ✅ Monitor PR for feedback
2. ✅ Respond to review comments
3. ✅ Make requested changes
4. ✅ Keep PR updated with upstream changes
5. ✅ Continue with implementation phases

---

## Quick Commands Summary

```bash
# Commit setup files
git add .vscode/ tools/scripts/audit_qt_dependencies.sh
git commit -m "Setup: Add VS Code configuration for Qt removal (Issue #8200)"

# Push to your fork
git push origin feature/remove-qt-from-core

# After that, create PR on GitHub web interface
# Link: https://github.com/Prachi194agrawal/OpenMS
```

Good luck with your PR! 🚀
