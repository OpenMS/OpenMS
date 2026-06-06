#!/usr/bin/env python3
"""Print imzML dataset and per-spectrum metadata after ImzMLFile.load().

Demonstrates that ImzMLMeta fields are mirrored as MetaValues on MSExperiment
and that pixel coordinates are attached to each MSSpectrum.

Run from the repo root (recommended):

    OpenMS-build/pyOpenMS/.venv/bin/python src/pyOpenMS/demo_imzML_metadata.py

Or pass one or more .imzML paths (companion .ibd must sit beside each file):

    OpenMS-build/pyOpenMS/.venv/bin/python src/pyOpenMS/demo_imzML_metadata.py \\
        path/to/file.imzML

Set PYOPENMS_BUILD if bindings live outside OpenMS-build/pyOpenMS.
"""

from __future__ import annotations

import argparse
import os
import sys


def _setup_pyopenms_path() -> None:
    src_pyopenms = os.path.dirname(os.path.abspath(__file__))
    if src_pyopenms in sys.path:
        sys.path.remove(src_pyopenms)

    build_pyopenms = os.environ.get(
        "PYOPENMS_BUILD",
        os.path.abspath(os.path.join(src_pyopenms, "..", "..", "OpenMS-build", "pyOpenMS")),
    )
    if os.path.isdir(os.path.join(build_pyopenms, "pyopenms")):
        if build_pyopenms not in sys.path:
            sys.path.insert(0, build_pyopenms)


_setup_pyopenms_path()

import pyopenms  # noqa: E402


DATASET_META_KEYS = (
    "imzml:imaging_mode",
    "imzml:max_count_x",
    "imzml:max_count_y",
    "imzml:max_count_z",
    "imzml:pixel_size_x",
    "imzml:pixel_size_y",
    "imzml:max_dim_x",
    "imzml:max_dim_y",
    "imzml:uuid",
    "imzml:ibd_md5",
    "imzml:ibd_sha1",
    "imzml:mz_data_type",
    "imzml:int_data_type",
    "imzml:scan_pattern",
    "imzml:scan_direction",
    "imzml:line_scan_direction",
    "imzml:polarity",
    "imzml:ibd_path",
)


def _repo_test_data_dir() -> str:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    return os.path.join(repo_root, "src", "tests", "class_tests", "openms", "data")


def _default_fixture_paths() -> list[str]:
    data = _repo_test_data_dir()
    return [
        os.path.join(data, "ImzMLFile_1_Example_Continuous.imzML"),
        os.path.join(data, "ImzMLFile_2_Example_Processed.imzML"),
    ]


def _infer_ibd_path(imzml_path: str) -> str:
    lower = imzml_path.lower()
    if lower.endswith(".imzml"):
        return imzml_path[:-6] + ".ibd"
    return imzml_path + ".ibd"


def _present(value: str) -> str:
    return value if value else "(not present)"


def print_dataset_metadata(label: str, imzml_path: str) -> None:
    print("=" * 72)
    print(f"  {label}")
    print(f"  File: {imzml_path}")
    print("=" * 72)

    exp = pyopenms.MSExperiment()
    pyopenms.ImzMLFile().load(imzml_path, exp)

    print("\n--- MSExperiment MetaValues (after ImzMLFile.load) ---")
    for key in DATASET_META_KEYS:
        if exp.metaValueExists(key):
            val = str(exp.getMetaValue(key))
            if key == "imzml:ibd_path" and len(val) > 72:
                val = "..." + val[-69:]
            print(f"  {key:30s}  {val}")
        else:
            print(f"  {key:30s}  (not present)")

    meta, index = pyopenms.ImzMLFile().loadSpectraIndex(imzml_path)
    print("\n--- ImzMLMeta struct (loadSpectraIndex) ---")
    print(f"  imaging_mode         {meta.imaging_mode}")
    print(f"  max_count_x/y/z      {meta.max_count_x} / {meta.max_count_y} / {meta.max_count_z}")
    print(f"  pixel_size_x/y (um)  {meta.pixel_size_x} / {meta.pixel_size_y}")
    print(f"  max_dim_x/y (um)     {meta.max_dim_x} / {meta.max_dim_y}")
    print(f"  uuid                 {_present(meta.uuid)}")
    print(f"  ibd_md5              {_present(meta.ibd_md5)}")
    print(f"  ibd_sha1             {_present(meta.ibd_sha1)}")
    print(f"  scan_pattern         {_present(meta.scan_pattern)}")
    print(f"  scan_direction       {_present(meta.scan_direction)}")
    print(f"  line_scan_direction  {_present(meta.line_scan_direction)}")
    print(f"  polarity             {_present(meta.polarity)}")
    print(f"  index entries        {len(index)}")

    if exp.size() == 0:
        print("\n  (no spectra loaded)")
        print()
        return

    spec = exp.getSpectrum(0)
    print("\n--- Per-spectrum MetaValues (spectrum[0]) ---")
    for key in ("imzml:x", "imzml:y", "imzml:z"):
        if spec.metaValueExists(key):
            print(f"  {key:20s}  {spec.getMetaValue(key)}")
    print(f"  peaks                {spec.size()}")
    if spec.size() > 0:
        mzs, ints = spec.get_peaks()
        print(f"  first peak m/z       {float(mzs[0]):.6f}")
        print(f"  first peak intensity {float(ints[0]):.6f}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print imzML dataset metadata mirrored on MSExperiment.",
    )
    parser.add_argument(
        "imzml",
        nargs="*",
        help="Path(s) to .imzML file(s); default: both reference fixtures",
    )
    args = parser.parse_args()

    paths = [os.path.abspath(p) for p in args.imzml] if args.imzml else _default_fixture_paths()

    for path in paths:
        if not os.path.isfile(path):
            print(f"File not found: {path}", file=sys.stderr)
            sys.exit(1)
        ibd = _infer_ibd_path(path)
        if not os.path.isfile(ibd):
            print(f"Companion .ibd not found: {ibd}", file=sys.stderr)
            sys.exit(1)

    labels = {
        "ImzMLFile_1_Example_Continuous.imzML": "CONTINUOUS imzML (IMS:1000030)",
        "ImzMLFile_2_Example_Processed.imzML": "PROCESSED imzML (IMS:1000031)",
    }

    for path in paths:
        name = os.path.basename(path)
        label = labels.get(name, f"imzML dataset ({name})")
        print_dataset_metadata(label, path)

    print("OpenMS imzML — dataset metadata on MSExperiment / MSSpectrum")


if __name__ == "__main__":
    main()
