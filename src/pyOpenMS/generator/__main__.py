"""
pyOpenMS Generator - Command Line Interface

This module provides the main entry point for the code generator.
Uses libclang for accurate C++ type information.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

from .pxd_parser import PxdParser
from .nanobind_emitter import NanobindEmitter, DOMAIN_NAMES
from .type_registry import TypeRegistry

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
        description="Generate nanobind C++ bindings from pyOpenMS .pxd files using libclang",
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
        required=False,
        default=None,
        help="(Deprecated) Directory containing .pyx addon files. No longer used.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for generated C++ files",
    )

    parser.add_argument(
        "--openms-include-dir",
        type=Path,
        nargs="+",
        required=True,
        help="OpenMS include directories for libclang parsing",
    )

    parser.add_argument(
        "--num-modules",
        type=int,
        default=None,
        help="Deprecated, ignored. Modules are now split by domain.",
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
        "--libclang-cache-dir",
        type=Path,
        help="Directory to cache libclang parse results (speeds up subsequent runs)",
    )

    parser.add_argument(
        "--libclang-batch-mode",
        choices=["auto", "on", "off"],
        default="auto",
        help="Batch-parse multiple headers per translation unit to avoid repeatedly parsing common includes",
    )

    parser.add_argument(
        "--libclang-batch-size",
        type=int,
        default=50,
        help="Number of headers per batch when libclang batching is enabled",
    )

    parser.add_argument(
        "--core-only",
        action="store_true",
        default=False,
        help="Deprecated, ignored. All classes are always bound.",
    )

    parser.add_argument(
        "--all-classes",
        action="store_true",
        default=False,
        help="Deprecated, ignored. All classes are always bound.",
    )

    return parser.parse_args(args)


def main(args: Optional[List[str]] = None) -> int:
    """Main entry point for the generator."""
    opts = parse_args(args)

    if opts.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("pyOpenMS Generator (libclang)")
    logger.info(f"  PXD directory: {opts.pxd_dir}")
    logger.info(f"  Output directory: {opts.output_dir}")
    logger.info(f"  Include directories: {opts.openms_include_dir}")
    logger.info(f"  Module split: domain-based ({len(DOMAIN_NAMES)} domains)")
    if opts.libclang_cache_dir:
        logger.info(f"  Cache directory: {opts.libclang_cache_dir}")
    logger.info(f"  Libclang batch mode: {opts.libclang_batch_mode} (size={opts.libclang_batch_size})")

    # Check libclang availability
    if not CLANG_AVAILABLE:
        logger.error("libclang is not available")
        logger.error("Install with: pip install clang==14.0")
        return 1

    # Validate input directories
    if not opts.pxd_dir.exists():
        logger.error(f"PXD directory does not exist: {opts.pxd_dir}")
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

    logger.info(f"Parsed {len(classes)} classes from .pxd files")

    # Note: Legacy .pyx addon processing was removed. Addons are now pure Python
    # in pyopenms/addons/ and loaded at runtime. C++ special methods are hardcoded
    # in nanobind_emitter.py's SPECIAL_METHODS dict.

    if opts.dry_run:
        logger.info("Dry run - not generating output")
        return 0

    # Parse C++ headers with libclang
    logger.info("Parsing C++ headers with libclang...")
    cpp_parser = CppHeaderParser(
        include_paths=opts.openms_include_dir,
        cache_dir=opts.libclang_cache_dir,
    )

    # Deduplicate headers (multiple classes can share a header)
    header_rels = []
    seen_headers = set()
    for class_decl in classes.values():
        if not class_decl.header_file:
            continue
        header_rel = class_decl.header_file.strip("<>")
        if header_rel in seen_headers:
            continue
        header_rels.append(header_rel)
        seen_headers.add(header_rel)

    header_paths = []
    missing_headers = []
    for header_rel in header_rels:
        found = False
        for inc_dir in opts.openms_include_dir:
            header_path = inc_dir / header_rel
            if header_path.exists():
                header_paths.append(header_path)
                found = True
                break
        if not found:
            missing_headers.append(header_rel)

    if missing_headers:
        logger.warning(f"Missing {len(missing_headers)} headers referenced by .pxd files")
        if opts.verbose:
            for h in missing_headers[:20]:
                logger.warning(f"  Missing header: {h}")

    cpp_classes = cpp_parser.parse_headers(
        header_paths,
        batch_mode=opts.libclang_batch_mode,
        batch_size=opts.libclang_batch_size,
    )

    logger.info(f"Parsed {len(cpp_classes)} C++ classes from headers")

    # Get scoped enums (C++11 enum class) from C++ headers
    scoped_enums = cpp_parser.get_scoped_enums()
    logger.info(f"Found {len(scoped_enums)} scoped enums (enum class)")

    # Merge C++ info with .pxd allowlist
    logger.info("Merging C++ and .pxd information...")
    merged_classes = merge_with_pxd(cpp_classes, classes, scoped_enums)
    logger.info(f"Merged {len(merged_classes)} classes")

    # Generate bindings
    logger.info("Generating nanobind C++ bindings...")
    if opts.num_modules is not None:
        logger.warning("--num-modules is deprecated and ignored. Using domain-based split.")
    emitter = NanobindEmitter()

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
