#!/usr/bin/env python3
"""
Automated script to migrate plain C++ enums to enum class in OpenMS.
This script handles:
1. Updating C++ headers (enum -> enum class)
2. Updating C++ source files (array sizes, loops)
3. Updating all usage sites in code
4. Updating Cython .pxd files
"""

import os
import re
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Set

# Base paths
OPENMS_ROOT = Path(__file__).parent
INCLUDE_DIR = OPENMS_ROOT / "src/openms/include/OpenMS"
SOURCE_DIR = OPENMS_ROOT / "src/openms/source"
PXD_DIR = OPENMS_ROOT / "src/pyOpenMS/pxds"

# Priority list of enums to migrate (from your analysis)
PRIORITY_ENUMS = [
    ("SpectrumSettings", "SpectrumType", "METADATA/SpectrumSettings.h"),
    ("ChromatogramSettings", "ChromatogramType", "METADATA/ChromatogramSettings.h"),
    ("IonSource", "Polarity", "METADATA/IonSource.h"),
    ("MassAnalyzer", "AnalyzerType", "METADATA/MassAnalyzer.h"),
    ("IonSource", "IonizationMethod", "METADATA/IonSource.h"),
    ("InstrumentSettings", "ScanMode", "METADATA/InstrumentSettings.h"),
    ("MassAnalyzer", "ResolutionMethod", "METADATA/MassAnalyzer.h"),
    ("IonDetector", "Type", "METADATA/IonDetector.h"),
    ("BaseFeature", "AnnotationState", "KERNEL/BaseFeature.h"),
    ("MassAnalyzer", "ScanDirection", "METADATA/MassAnalyzer.h"),
    ("MassAnalyzer", "ResolutionType", "METADATA/MassAnalyzer.h"),
    ("MassAnalyzer", "ScanLaw", "METADATA/MassAnalyzer.h"),
    ("MassAnalyzer", "ReflectronState", "METADATA/MassAnalyzer.h"),
    ("IonSource", "InletType", "METADATA/IonSource.h"),
    ("IonDetector", "AcquisitionMode", "METADATA/IonDetector.h"),
    ("ProteinIdentification", "PeakMassType", "METADATA/ProteinIdentification.h"),
    ("SourceFile", "ChecksumType", "METADATA/SourceFile.h"),
    ("Instrument", "IonOpticsType", "METADATA/Instrument.h"),
]


class EnumInfo:
    """Holds information about an enum to migrate"""
    def __init__(self, class_name: str, enum_name: str, header_path: str):
        self.class_name = class_name
        self.enum_name = enum_name
        self.header_path = header_path
        self.enum_values: List[str] = []
        self.pxd_path: str = ""
        self.cpp_path: str = ""

    def __repr__(self):
        return f"{self.class_name}::{self.enum_name}"


def find_enum_values(header_file: Path, enum_name: str) -> List[str]:
    """Extract enum value names from header file"""
    values = []

    if not header_file.exists():
        return values

    with open(header_file, 'r') as f:
        content = f.read()

    # Find the enum block - handle both regular and inline declarations with newlines
    pattern = r'enum\s+(?:class\s+)?' + enum_name + r'\s*\n?\s*\{([^}]+)\}'
    match = re.search(pattern, content, re.DOTALL)

    if match:
        enum_body = match.group(1)
        # Extract individual values (handle comments and assignments)
        for line in enum_body.split(','):
            # Remove comments and whitespace
            line = re.sub(r'//.*', '', line)
            line = re.sub(r'/\*.*?\*/', '', line, flags=re.DOTALL)
            line = line.strip()

            if line:
                # Get just the identifier (before any = assignment)
                value = re.match(r'([A-Z_][A-Z0-9_]*)', line)
                if value:
                    values.append(value.group(1))

    return values


def update_header_file(enum_info: EnumInfo) -> bool:
    """Update C++ header to use enum class"""
    # header_path is like "METADATA/SpectrumSettings.h" - use it directly
    header_file = INCLUDE_DIR / enum_info.header_path

    if not header_file.exists():
        print(f"  ❌ Header not found: {header_file}")
        return False

    with open(header_file, 'r') as f:
        content = f.read()

    original_content = content

    # Pattern 1: enum EnumName -> enum class EnumName (handle newlines)
    pattern = r'\benum\s+' + enum_info.enum_name + r'\s*\n?\s*\{'
    replacement = f'enum class {enum_info.enum_name}\n    {{'
    content = re.sub(pattern, replacement, content)

    # Pattern 2: Update array size declarations
    # NamesOfEnumName[SIZE_OF_ENUMNAME] -> NamesOfEnumName[static_cast<size_t>(EnumName::SIZE_OF_ENUMNAME)]
    if enum_info.enum_values:
        size_constant = enum_info.enum_values[-1]  # Usually SIZE_OF_XXX
        if size_constant.startswith("SIZE_OF_"):
            pattern = rf'(\w+)\[{size_constant}\]'
            replacement = rf'\1[static_cast<size_t>({enum_info.enum_name}::{size_constant})]'
            content = re.sub(pattern, replacement, content)

    if content != original_content:
        with open(header_file, 'w') as f:
            f.write(content)
        print(f"  ✅ Updated header: {header_file.name}")
        return True
    else:
        print(f"  ⚠️  No changes needed in header: {header_file.name}")
        return False


def update_source_file(enum_info: EnumInfo) -> bool:
    """Update C++ source file"""
    # Try to find corresponding .cpp file - use header_path directly
    cpp_path = enum_info.header_path.replace(".h", ".cpp")
    source_file = SOURCE_DIR / cpp_path

    if not source_file.exists():
        print(f"  ℹ️  No source file found: {source_file.name}")
        return False

    with open(source_file, 'r') as f:
        content = f.read()

    original_content = content

    # Update SIZE_OF_XXX references
    if enum_info.enum_values:
        size_constant = enum_info.enum_values[-1]
        if size_constant.startswith("SIZE_OF_"):
            # Pattern 1: for loops
            pattern = rf'\b{size_constant}\b'
            replacement = f'static_cast<size_t>({enum_info.enum_name}::{size_constant})'
            content = re.sub(pattern, replacement, content)

            # Pattern 2: Constructor initializers (unqualified enum values)
            for value in enum_info.enum_values[:-1]:  # Exclude SIZE_OF_
                # Look for unqualified usage in constructors
                pattern = rf'\b({value})\b(?!\s*::|,)'
                # Check if it's not already qualified
                if re.search(rf'(?<!::){value}\b(?!::)', content):
                    replacement = f'{enum_info.enum_name}::{value}'
                    content = re.sub(pattern, replacement, content)

    if content != original_content:
        with open(source_file, 'w') as f:
            f.write(content)
        print(f"  ✅ Updated source: {source_file.name}")
        return True
    else:
        print(f"  ⚠️  No changes needed in source: {source_file.name}")
        return False


def update_usage_sites(enum_info: EnumInfo) -> int:
    """Update all usage sites in C++ code"""
    count = 0

    # Find all C++ and header files
    for ext in ['*.cpp', '*.h']:
        for filepath in SOURCE_DIR.rglob(ext):
            with open(filepath, 'r') as f:
                content = f.read()

            original_content = content

            # Update each enum value: ClassName::VALUE -> ClassName::EnumName::VALUE
            for value in enum_info.enum_values:
                # Skip SIZE_OF_ constants for now
                if value.startswith("SIZE_OF_"):
                    continue

                # Pattern: ClassName::VALUE but not ClassName::EnumName::VALUE
                pattern = rf'\b{enum_info.class_name}::{value}\b(?!::)'
                replacement = f'{enum_info.class_name}::{enum_info.enum_name}::{value}'
                content = re.sub(pattern, replacement, content)

            if content != original_content:
                with open(filepath, 'w') as f:
                    f.write(content)
                count += 1

    # Also check test files
    test_dir = OPENMS_ROOT / "src/tests/class_tests/openms/source"
    if test_dir.exists():
        for filepath in test_dir.glob("*.cpp"):
            with open(filepath, 'r') as f:
                content = f.read()

            original_content = content

            for value in enum_info.enum_values:
                if value.startswith("SIZE_OF_"):
                    continue
                pattern = rf'\b{enum_info.class_name}::{value}\b(?!::)'
                replacement = f'{enum_info.class_name}::{enum_info.enum_name}::{value}'
                content = re.sub(pattern, replacement, content)

            if content != original_content:
                with open(filepath, 'w') as f:
                    f.write(content)
                count += 1

    if count > 0:
        print(f"  ✅ Updated {count} usage sites")
    else:
        print(f"  ℹ️  No usage sites needed updating")

    return count


def update_pxd_file(enum_info: EnumInfo) -> bool:
    """Update Cython .pxd file"""
    # Find .pxd file
    pxd_file = PXD_DIR / f"{enum_info.class_name}.pxd"

    if not pxd_file.exists():
        print(f"  ❌ PXD file not found: {pxd_file}")
        return False

    with open(pxd_file, 'r') as f:
        content = f.read()

    original_content = content

    # Pattern: cdef enum EnumName: -> cdef enum class EnumName "OpenMS::ClassName::EnumName":
    pattern = rf'cdef enum {enum_info.enum_name}:'
    replacement = f'cdef enum class {enum_info.enum_name} "OpenMS::{enum_info.class_name}::{enum_info.enum_name}":'
    content = re.sub(pattern, replacement, content)

    if content != original_content:
        with open(pxd_file, 'w') as f:
            f.write(content)
        print(f"  ✅ Updated PXD: {pxd_file.name}")
        return True
    else:
        print(f"  ⚠️  No changes needed in PXD: {pxd_file.name}")
        return False


def migrate_enum(class_name: str, enum_name: str, header_path: str) -> bool:
    """Migrate a single enum to enum class"""
    print(f"\n{'='*60}")
    print(f"Migrating: {class_name}::{enum_name}")
    print(f"{'='*60}")

    enum_info = EnumInfo(class_name, enum_name, header_path)

    # Find enum values - header_path is like "METADATA/SpectrumSettings.h"
    header_file = INCLUDE_DIR / header_path
    enum_info.enum_values = find_enum_values(header_file, enum_name)

    if not enum_info.enum_values:
        print(f"  ❌ Could not find enum values in {header_file}")
        return False

    print(f"  Found {len(enum_info.enum_values)} enum values: {', '.join(enum_info.enum_values[:5])}...")

    # Step 1: Update header
    update_header_file(enum_info)

    # Step 2: Update source file
    update_source_file(enum_info)

    # Step 3: Update usage sites
    update_usage_sites(enum_info)

    # Step 4: Update .pxd file
    update_pxd_file(enum_info)

    print(f"  ✅ Migration complete for {class_name}::{enum_name}")
    return True


def main():
    """Main migration workflow"""
    print("=" * 60)
    print("OpenMS Enum to Enum Class Migration Script")
    print("=" * 60)
    print(f"Root directory: {OPENMS_ROOT}")
    print(f"Total enums to migrate: {len(PRIORITY_ENUMS)}")

    success_count = 0
    failed_enums = []

    for class_name, enum_name, header_path in PRIORITY_ENUMS:
        try:
            if migrate_enum(class_name, enum_name, header_path):
                success_count += 1
            else:
                failed_enums.append(f"{class_name}::{enum_name}")
        except Exception as e:
            print(f"  ❌ Error migrating {class_name}::{enum_name}: {e}")
            failed_enums.append(f"{class_name}::{enum_name}")

    # Summary
    print("\n" + "=" * 60)
    print("MIGRATION SUMMARY")
    print("=" * 60)
    print(f"✅ Successfully migrated: {success_count}/{len(PRIORITY_ENUMS)}")

    if failed_enums:
        print(f"\n❌ Failed migrations:")
        for enum in failed_enums:
            print(f"  - {enum}")

    print("\n" + "=" * 60)
    print("NEXT STEPS:")
    print("=" * 60)
    print("1. Review the changes with git diff")
    print("2. Build C++ code: cmake --build build")
    print("3. Run tests: ctest")
    print("4. Build pyOpenMS: cd src/pyOpenMS && python create_cpp_extension.py build")
    print("5. Run Python tests: pytest src/pyOpenMS/tests")


if __name__ == "__main__":
    main()
