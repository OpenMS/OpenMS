"""
pyOpenMS2 Generator - Command Line Interface

This module provides the main entry point for the code generator.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

from .pxd_parser import PxdParser
from .nanobind_emitter_v2 import NanobindEmitterV2
from .type_registry import TypeRegistry
from .addon_processor import AddonProcessor
from .cpp_parser import pxd_to_merged
from .doxygen_parser import DoxygenXMLParser, merge_doxygen_with_pxd

# Check if libclang is available
try:
    from .cpp_parser import CLANG_AVAILABLE, CppHeaderParser, merge_with_pxd
except ImportError:
    CLANG_AVAILABLE = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate nanobind C++ bindings from pyOpenMS .pxd files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--pxd-dir",
        type=Path,
        required=True,
        help="Directory containing .pxd files",
    )

    parser.add_argument(
        "--addons-dir",
        type=Path,
        required=True,
        help="Directory containing .pyx addon files",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for generated C++ files",
    )

    parser.add_argument(
        "--num-modules",
        type=int,
        default=8,
        help="Number of modules to split bindings into",
    )

    parser.add_argument(
        "--classes",
        nargs="*",
        help="Only generate bindings for specific classes (for testing)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse files but don't generate output",
    )

    parser.add_argument(
        "--use-libclang",
        action="store_true",
        help="Use libclang to parse C++ headers for accurate type information",
    )

    parser.add_argument(
        "--openms-include-dir",
        type=Path,
        nargs="+",
        help="OpenMS include directories for libclang parsing",
    )

    parser.add_argument(
        "--libclang-cache-dir",
        type=Path,
        help="Directory to cache libclang parse results (speeds up subsequent runs)",
    )

    parser.add_argument(
        "--core-only",
        action="store_true",
        default=True,
        help="Only bind core classes (default: True)",
    )

    parser.add_argument(
        "--all-classes",
        action="store_true",
        help="Attempt to bind all classes (may have compile errors)",
    )

    parser.add_argument(
        "--use-doxygen",
        action="store_true",
        help="Use Doxygen XML for accurate C++ type information",
    )

    parser.add_argument(
        "--doxygen-xml-dir",
        type=Path,
        help="Directory containing Doxygen XML output (required with --use-doxygen)",
    )

    return parser.parse_args(args)


def main(args: Optional[List[str]] = None) -> int:
    """Main entry point for the generator."""
    opts = parse_args(args)

    if opts.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info(f"pyOpenMS2 Generator")
    logger.info(f"  PXD directory: {opts.pxd_dir}")
    logger.info(f"  Addons directory: {opts.addons_dir}")
    logger.info(f"  Output directory: {opts.output_dir}")
    logger.info(f"  Number of modules: {opts.num_modules}")

    # Validate input directories
    if not opts.pxd_dir.exists():
        logger.error(f"PXD directory does not exist: {opts.pxd_dir}")
        return 1

    if not opts.addons_dir.exists():
        logger.error(f"Addons directory does not exist: {opts.addons_dir}")
        return 1

    # Create output directory
    opts.output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize type registry
    logger.info("Initializing type registry...")
    type_registry = TypeRegistry()

    # Parse .pxd files
    logger.info("Parsing .pxd files...")
    parser = PxdParser(type_registry)

    pxd_files = sorted(opts.pxd_dir.glob("*.pxd"))
    logger.info(f"Found {len(pxd_files)} .pxd files")

    classes = {}
    for pxd_file in pxd_files:
        # Filter if specific classes requested
        if opts.classes:
            stem = pxd_file.stem
            if stem not in opts.classes:
                continue

        logger.debug(f"Parsing {pxd_file.name}...")
        try:
            parsed = parser.parse_file(pxd_file)
            if parsed:
                classes.update(parsed)
        except Exception as e:
            logger.warning(f"Failed to parse {pxd_file.name}: {e}")
            if opts.verbose:
                import traceback

                traceback.print_exc()

    logger.info(f"Parsed {len(classes)} classes")

    # Process addon files
    logger.info("Processing addon files...")
    addon_processor = AddonProcessor()

    addon_files = sorted(opts.addons_dir.glob("*.pyx"))
    logger.info(f"Found {len(addon_files)} addon files")

    addons = {}
    for addon_file in addon_files:
        logger.debug(f"Processing addon {addon_file.name}...")
        try:
            parsed = addon_processor.parse_file(addon_file)
            if parsed:
                class_name = addon_file.stem
                addons[class_name] = parsed
        except Exception as e:
            logger.warning(f"Failed to process addon {addon_file.name}: {e}")
            if opts.verbose:
                import traceback

                traceback.print_exc()

    logger.info(f"Processed {len(addons)} addon files")

    if opts.dry_run:
        logger.info("Dry run - not generating output")
        return 0

    # Check if Doxygen mode is requested
    if opts.use_doxygen:
        if not opts.doxygen_xml_dir:
            logger.error("--doxygen-xml-dir required when using --use-doxygen")
            return 1

        if not opts.doxygen_xml_dir.exists():
            logger.error(f"Doxygen XML directory does not exist: {opts.doxygen_xml_dir}")
            return 1

        logger.info("Using Doxygen XML for accurate C++ type information")
        logger.info(f"  Doxygen XML directory: {opts.doxygen_xml_dir}")

        # Parse Doxygen XML
        logger.info("Parsing Doxygen XML files...")
        doxy_parser = DoxygenXMLParser(opts.doxygen_xml_dir)
        doxy_classes = doxy_parser.parse_all()
        logger.info(f"Parsed {len(doxy_classes)} classes from Doxygen XML")

        # Merge with .pxd allowlist
        logger.info("Merging Doxygen info with .pxd allowlist...")
        merged_classes = merge_doxygen_with_pxd(doxy_classes, classes)
        logger.info(f"Merged {len(merged_classes)} classes")

    # Check if libclang mode is requested
    elif opts.use_libclang:
        if not CLANG_AVAILABLE:
            logger.error("--use-libclang requested but libclang is not available")
            logger.error("Install with: pip install clang==14.0")
            return 1

        if not opts.openms_include_dir:
            logger.error("--openms-include-dir required when using --use-libclang")
            return 1

        logger.info("Using libclang for accurate C++ type information")
        logger.info(f"  Include directories: {opts.openms_include_dir}")
        if opts.libclang_cache_dir:
            logger.info(f"  Cache directory: {opts.libclang_cache_dir}")

        # Parse C++ headers
        logger.info("Parsing C++ headers with libclang...")
        cpp_parser = CppHeaderParser(
            include_paths=opts.openms_include_dir,
            cache_dir=opts.libclang_cache_dir,
        )
        cpp_classes = {}

        # Build mapping from class name to header file
        for class_name, class_decl in classes.items():
            if class_decl.header_file:
                # Strip angle brackets from header path (e.g., "<OpenMS/KERNEL/Peak1D.h>")
                header_rel = class_decl.header_file.strip("<>")
                # Find the header file in include directories
                for inc_dir in opts.openms_include_dir:
                    header_path = inc_dir / header_rel
                    if header_path.exists():
                        logger.debug(f"Parsing header for {class_name}: {header_path}")
                        try:
                            parsed = cpp_parser.parse_header(header_path)
                            cpp_classes.update(parsed)
                        except Exception as e:
                            logger.warning(f"Failed to parse {header_path}: {e}")
                        break

        logger.info(f"Parsed {len(cpp_classes)} C++ classes from headers")

        # Merge C++ info with .pxd allowlist
        logger.info("Merging C++ and .pxd information...")
        merged_classes = merge_with_pxd(cpp_classes, classes)
        logger.info(f"Merged {len(merged_classes)} classes")
    else:
        # Convert .pxd classes to MergedClass format for v2 emitter
        logger.info("Converting .pxd classes to merged format...")
        merged_classes = pxd_to_merged(classes)
        logger.info(f"Converted {len(merged_classes)} classes")

    # Generate with v2 emitter (always used now)
    core_only = opts.core_only and not opts.all_classes
    logger.info(f"Generating nanobind C++ bindings (core_only={core_only})...")
    emitter = NanobindEmitterV2(num_modules=opts.num_modules, core_only=core_only)

    try:
        emitter.emit(merged_classes, opts.output_dir)
    except Exception as e:
        logger.error(f"Failed to generate bindings: {e}")
        if opts.verbose:
            import traceback
            traceback.print_exc()
        return 1

    logger.info("Generation complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
