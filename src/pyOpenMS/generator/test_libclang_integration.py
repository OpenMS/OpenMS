#!/usr/bin/env python3
"""
Test script for libclang integration.

This script:
1. Parses a C++ header with libclang
2. Parses corresponding .pxd file
3. Merges them
4. Generates a sample binding file
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from generator.cpp_parser import CLANG_AVAILABLE, CppHeaderParser, merge_with_pxd
from generator.pxd_parser import PxdParser
from generator.type_registry import TypeRegistry
from generator.nanobind_emitter import NanobindEmitter


def main():
    if not CLANG_AVAILABLE:
        print("Error: libclang not available. Install with: pip install clang==14.0")
        sys.exit(1)

    openms_src = Path("/home/sachsenb/Development/tmp/OpenMS/src")
    openms_build = Path("/home/sachsenb/Development/tmp/OpenMS-build/src")
    output_dir = Path("/tmp/pyopenms_test")

    # Include paths for C++ parsing
    include_paths = [
        openms_src / "openms/include",
        openms_build / "openms/include",
        Path("/usr/include/x86_64-linux-gnu/qt6"),
        Path("/usr/include/x86_64-linux-gnu/qt6/QtCore"),
    ]

    # Test classes to generate
    test_classes = [
        ("Peak1D", "OpenMS/KERNEL/Peak1D.h", "Peak1D.pxd"),
        ("Peak2D", "OpenMS/KERNEL/Peak2D.h", "Peak2D.pxd"),
        ("ChromatogramPeak", "OpenMS/KERNEL/ChromatogramPeak.h", "ChromatogramPeak.pxd"),
    ]

    print("=== Parsing C++ headers with libclang ===")
    cpp_parser = CppHeaderParser(include_paths=include_paths)
    cpp_classes = {}

    for class_name, header_path, _ in test_classes:
        full_path = openms_src / "openms/include" / header_path
        if not full_path.exists():
            print(f"  Warning: Header not found: {full_path}")
            continue
        print(f"  Parsing {header_path}...")
        parsed = cpp_parser.parse_header(full_path)
        cpp_classes.update(parsed)

    print(f"  Found {len(cpp_classes)} C++ classes")
    for name, cls in cpp_classes.items():
        if "OpenMS" in cls.qualified_name:
            print(f"    - {cls.qualified_name}: {len(cls.methods)} methods, {len(cls.constructors)} ctors")
            print(f"      const methods: {sum(1 for m in cls.methods if m.is_const)}")
            print(f"      static methods: {sum(1 for m in cls.methods if m.is_static)}")

    print("\n=== Parsing .pxd files ===")
    type_registry = TypeRegistry()
    pxd_parser = PxdParser(type_registry)
    pxd_dir = openms_src / "pyOpenMS/pxds"
    pxd_classes = {}

    for class_name, _, pxd_file in test_classes:
        full_path = pxd_dir / pxd_file
        if not full_path.exists():
            print(f"  Warning: .pxd file not found: {full_path}")
            continue
        print(f"  Parsing {pxd_file}...")
        parsed = pxd_parser.parse_file(full_path)
        pxd_classes.update(parsed)

    print(f"  Found {len(pxd_classes)} .pxd classes")
    for name, cls in pxd_classes.items():
        print(f"    - {name}: {len(cls.methods)} methods")

    print("\n=== Merging C++ and .pxd info ===")
    merged = merge_with_pxd(cpp_classes, pxd_classes)

    print(f"  Merged {len(merged)} classes")
    for name, cls in merged.items():
        print(f"    - {name}:")
        print(f"      Methods to wrap: {len(cls.methods)}")
        print(f"      Constructors: {len(cls.constructors)}")
        for mm in cls.methods[:5]:
            m = mm.cpp_method
            const_str = " const" if m.is_const else ""
            params = ", ".join(p.type_str for p in m.parameters)
            print(f"        {m.return_type} {m.name}({params}){const_str}")
        if len(cls.methods) > 5:
            print(f"        ... and {len(cls.methods) - 5} more")

    print("\n=== Generating bindings ===")
    import shutil
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    emitter = NanobindEmitter()
    emitter.emit(merged, output_dir)

    print(f"  Generated files in {output_dir}:")
    for f in output_dir.iterdir():
        print(f"    - {f.name}")
        if f.suffix == ".cpp":
            content = f.read_text()
            lines = content.split("\n")
            print(f"      ({len(lines)} lines)")
            # Show a snippet
            for line in lines[20:40]:
                print(f"      {line}")

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
