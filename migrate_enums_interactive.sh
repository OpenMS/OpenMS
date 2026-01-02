#!/bin/bash
# Interactive enum migration script for OpenMS
# Migrates enums one at a time with verification

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}OpenMS Enum Migration Script${NC}"
echo -e "${BLUE}========================================${NC}"

# Run the Python migration script
echo -e "\n${YELLOW}Running automated migration...${NC}\n"
python3 migrate_enums.py

# Verify changes
echo -e "\n${YELLOW}Checking for compilation errors...${NC}"
if [ -d "build" ]; then
    echo -e "${GREEN}Building C++ code...${NC}"
    cmake --build build 2>&1 | tee build_output.log

    if grep -q "error:" build_output.log; then
        echo -e "${RED}❌ Compilation errors found!${NC}"
        echo "Check build_output.log for details"
        exit 1
    else
        echo -e "${GREEN}✅ C++ code compiles successfully!${NC}"
    fi
else
    echo -e "${YELLOW}⚠️  No build directory found. Skipping compilation check.${NC}"
    echo "Run: mkdir build && cd build && cmake .. to create build directory"
fi

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}Migration complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "\nNext steps:"
echo -e "1. Review changes: ${YELLOW}git diff${NC}"
echo -e "2. Run tests: ${YELLOW}cd build && ctest${NC}"
echo -e "3. Build pyOpenMS: ${YELLOW}cd src/pyOpenMS && python create_cpp_extension.py build${NC}"
echo -e "4. Commit changes: ${YELLOW}git add -A && git commit -m 'Migrate enums to enum class'${NC}"
