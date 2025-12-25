#!/bin/bash
# Script to audit Qt dependencies in OpenMS and OpenSwath libraries
# This helps track progress on Issue #8200

set -e

WORKSPACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$WORKSPACE_ROOT"

echo "================================================================"
echo "Qt Dependency Audit for OpenMS Core Libraries"
echo "Issue: #8200 - Remove all Qt dependencies from OpenMS and OpenSwath"
echo "================================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}[1] Searching for Qt includes in src/openms/${NC}"
echo "------------------------------------------------"
QT_INCLUDES=$(grep -r "#include.*<Q[A-Z]" --include="*.h" --include="*.cpp" src/openms/ 2>/dev/null | grep -v "src/openms/extern" | grep -v "src/openms_gui" || true)
if [ -z "$QT_INCLUDES" ]; then
  echo -e "${GREEN}✓ No Qt includes found!${NC}"
else
  echo -e "${RED}Found Qt includes:${NC}"
  echo "$QT_INCLUDES" | head -20
  COUNT=$(echo "$QT_INCLUDES" | wc -l)
  echo ""
  echo -e "${RED}Total: $COUNT files with Qt includes${NC}"
fi
echo ""

echo -e "${YELLOW}[2] Searching for Qt includes in src/openswathalgo/${NC}"
echo "------------------------------------------------"
QT_INCLUDES_SWATH=$(grep -r "#include.*<Q[A-Z]" --include="*.h" --include="*.cpp" src/openswathalgo/ 2>/dev/null || true)
if [ -z "$QT_INCLUDES_SWATH" ]; then
  echo -e "${GREEN}✓ No Qt includes found in OpenSwathAlgo!${NC}"
else
  echo -e "${RED}Found Qt includes:${NC}"
  echo "$QT_INCLUDES_SWATH"
fi
echo ""

echo -e "${YELLOW}[3] Analyzing Qt class usage in OpenMS core${NC}"
echo "------------------------------------------------"
declare -A qt_classes=(
  ["QString"]="QString"
  ["QStringList"]="QStringList"
  ["QFile"]="QFile"
  ["QFileInfo"]="QFileInfo"
  ["QDir"]="QDir"
  ["QDateTime"]="QDateTime"
  ["QByteArray"]="QByteArray"
  ["QList"]="QList"
  ["QVector"]="QVector"
  ["QMap"]="QMap"
  ["QNetworkAccessManager"]="QNetworkAccessManager"
  ["QTcpSocket"]="QTcpSocket"
  ["QProcess"]="QProcess"
  ["QFileSystemWatcher"]="QFileSystemWatcher"
  ["QUrl"]="QUrl"
)

for class in "${!qt_classes[@]}"; do
  COUNT=$(grep -r "\b$class\b" --include="*.h" --include="*.cpp" src/openms/ 2>/dev/null | grep -v "src/openms/extern" | grep -v "src/openms_gui" | wc -l || echo "0")
  if [ "$COUNT" -gt 0 ]; then
    echo -e "${RED}  $class: $COUNT occurrences${NC}"
  fi
done
echo ""

echo -e "${YELLOW}[4] Files with Qt dependencies (excluding GUI)${NC}"
echo "------------------------------------------------"
FILES_WITH_QT=$(grep -rl "QString\|QFile\|QDir\|QNetwork\|QProcess" --include="*.h" --include="*.cpp" src/openms/source/ src/openms/include/ 2>/dev/null | grep -v "src/openms_gui" | sort -u || true)
if [ -z "$FILES_WITH_QT" ]; then
  echo -e "${GREEN}✓ No files with Qt dependencies found!${NC}"
else
  echo -e "${RED}Files needing refactoring:${NC}"
  echo "$FILES_WITH_QT" | nl
  FILE_COUNT=$(echo "$FILES_WITH_QT" | wc -l)
  echo ""
  echo -e "${RED}Total: $FILE_COUNT files need Qt removal${NC}"
fi
echo ""

echo -e "${YELLOW}[5] CMake Qt dependencies${NC}"
echo "------------------------------------------------"
CMAKE_QT=$(grep -r "Qt[56]::\|find_package.*Qt\|target_link_libraries.*Qt" --include="CMakeLists.txt" src/openms/ 2>/dev/null || true)
if [ -z "$CMAKE_QT" ]; then
  echo -e "${GREEN}✓ No Qt CMake dependencies found!${NC}"
else
  echo -e "${RED}CMake files with Qt:${NC}"
  echo "$CMAKE_QT"
fi
echo ""

echo -e "${YELLOW}[6] Suggested Replacements${NC}"
echo "------------------------------------------------"
cat << EOF
  QString          → std::string
  QStringList      → std::vector<std::string>
  QFile/QFileInfo  → std::filesystem::path (C++17)
  QDir             → std::filesystem::directory_iterator
  QDateTime        → std::chrono or boost::date_time
  QByteArray       → std::vector<char> or std::string
  QList/QVector    → std::vector
  QMap             → std::map or std::unordered_map
  QProcess         → custom wrapper using fork/exec or boost::process
  QNetworkAccessManager → boost::asio or cpp-httplib (see PR #8201)
  QFileSystemWatcher → inotify (Linux) or boost::asio::file_descriptor
EOF
echo ""

echo "================================================================"
echo -e "${GREEN}Audit Complete!${NC}"
echo "================================================================"
echo ""
echo "Next steps:"
echo "1. Create a feature branch: git checkout -b feature/remove-qt-from-core"
echo "2. Start with simple replacements (QString → std::string)"
echo "3. Test after each major change"
echo "4. Update CMakeLists.txt to remove Qt dependencies"
echo ""
