#!/bin/bash
# Quick setup script for working on Issue #8200
# Run this after opening the workspace in VS Code

set -e

WORKSPACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$WORKSPACE_ROOT"

echo "============================================"
echo "OpenMS Qt Removal Setup (Issue #8200)"
echo "============================================"
echo ""

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ] || [ ! -d "src/openms" ]; then
    echo "❌ Error: Not in OpenMS root directory"
    echo "Please run this script from the OpenMS workspace root"
    exit 1
fi

echo "✅ Found OpenMS workspace"
echo ""

# Check git status
echo "📊 Checking git status..."
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
echo "   Current branch: $CURRENT_BRANCH"
echo ""

# Check if develop branch exists
if git rev-parse --verify develop >/dev/null 2>&1; then
    echo "✅ Found develop branch"
else
    echo "⚠️  Warning: No 'develop' branch found"
    echo "   You may need to: git fetch upstream develop"
fi
echo ""

# Prompt to create feature branch
if [ "$CURRENT_BRANCH" != "feature/remove-qt-from-core" ]; then
    echo "📝 Current branch is not 'feature/remove-qt-from-core'"
    read -p "Create feature branch now? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git checkout -b feature/remove-qt-from-core 2>/dev/null || echo "Branch may already exist"
        echo "✅ Switched to feature/remove-qt-from-core"
    fi
fi
echo ""

# Run the Qt audit
echo "🔍 Running Qt dependency audit..."
echo ""
if [ -f "tools/scripts/audit_qt_dependencies.sh" ]; then
    bash tools/scripts/audit_qt_dependencies.sh
else
    echo "❌ Audit script not found at tools/scripts/audit_qt_dependencies.sh"
fi
echo ""

# Check for build directory
echo "🔨 Checking build configuration..."
if [ -d "build" ]; then
    echo "✅ Build directory exists"
    if [ -f "build/CMakeCache.txt" ]; then
        WITH_GUI=$(grep "WITH_GUI:BOOL" build/CMakeCache.txt | cut -d= -f2 || echo "UNKNOWN")
        echo "   WITH_GUI=$WITH_GUI"
    fi
else
    echo "📦 No build directory found"
    read -p "Configure build now? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Configuring with -DWITH_GUI=OFF for faster development builds..."
        cmake -B build -S . -DWITH_GUI=OFF -DCMAKE_BUILD_TYPE=Debug
        echo "✅ Build configured"
    fi
fi
echo ""

# Summary and next steps
echo "============================================"
echo "📚 Resources Available"
echo "============================================"
echo ""
echo "In .vscode/ directory:"
echo "  • README.md - Getting started guide"
echo "  • QT_REMOVAL_ACTION_PLAN.md - Detailed implementation plan"
echo "  • QT_REMOVAL_CHECKLIST.md - Track your progress"
echo "  • QT_REPLACEMENT_PATTERNS.md - Code examples"
echo "  • ARCHITECTURE_DIAGRAM.md - Visual overview"
echo "  • GITHUB_ISSUE_COMMENT.md - Post this on GitHub"
echo ""
echo "============================================"
echo "🚀 Next Steps"
echo "============================================"
echo ""
echo "1. Read the action plan:"
echo "   less .vscode/QT_REMOVAL_ACTION_PLAN.md"
echo ""
echo "2. Post on GitHub Issue #8200:"
echo "   Copy content from .vscode/GITHUB_ISSUE_COMMENT.md"
echo "   https://github.com/OpenMS/OpenMS/issues/8200"
echo ""
echo "3. Start with Phase 1 (String & DataValue):"
echo "   vim src/openms/include/OpenMS/DATASTRUCTURES/String.h"
echo ""
echo "4. Build and test:"
echo "   cmake --build build --target OpenMS -j4"
echo "   ctest --test-dir build -R String_test"
echo ""
echo "5. Track progress with checklist:"
echo "   code .vscode/QT_REMOVAL_CHECKLIST.md"
echo ""
echo "============================================"
echo "💡 Helpful Commands"
echo "============================================"
echo ""
echo "# Re-run audit"
echo "bash tools/scripts/audit_qt_dependencies.sh"
echo ""
echo "# Build OpenMS core"
echo "cmake --build build --target OpenMS -j4"
echo ""
echo "# Run specific test"
echo "ctest --test-dir build -R String_test -V"
echo ""
echo "# Check Qt linkage (should be empty after refactoring)"
echo "ldd build/lib/libOpenMS.so 2>/dev/null | grep -i qt || echo 'No Qt found ✅'"
echo ""
echo "============================================"
echo "Setup complete! Happy coding! 🎉"
echo "============================================"
