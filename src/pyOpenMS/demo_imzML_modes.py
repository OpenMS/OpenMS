#!/usr/bin/env python3
"""Demonstrate all imzML load modes in pyOpenMS (mirrors Tutorial_ImzMLFile.cpp).

Run from the pyOpenMS build directory (recommended):

    cd OpenMS-build/pyOpenMS
    ./.venv/bin/python ../../src/pyOpenMS/demo_imzML_modes.py \\
        ../../src/tests/class_tests/openms/data/ImzMLFile_1_Example_Continuous.imzML

Or pass any path to a .imzML file (companion .ibd must sit beside it).
"""

from __future__ import annotations

import argparse
import os
import sys


def _setup_pyopenms_path() -> None:
    """Prefer pyOpenMS bindings from OpenMS-build over the source tree."""
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

import pyopenms


def _default_continuous_path() -> str:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    return os.path.join(
        repo_root,
        "src",
        "tests",
        "class_tests",
        "openms",
        "data",
        "ImzMLFile_1_Example_Continuous.imzML",
    )


def print_peaks(spec: pyopenms.MSSpectrum, max_peaks: int = 3) -> None:
    n = min(spec.size(), max_peaks)
    for i in range(n):
        peak = spec[i]
        print(f"    m/z={peak.getMZ()}  I={peak.getIntensity()}")
    if spec.size() > n:
        print(f"    ... ({spec.size()} peaks total)")


def _infer_ibd_path(imzml_path: str) -> str:
    lower = imzml_path.lower()
    if lower.endswith(".imzml"):
        return imzml_path[:-6] + ".ibd"
    return imzml_path + ".ibd"


def _imzml_coords_from_index(od: pyopenms.OnDiscImzMLExperiment, index: int) -> tuple[int, int, int]:
    entry = od.getIndex(index)
    return int(entry.x), int(entry.y), int(entry.z)


def _geometry_coords(imzml_x: int, imzml_y: int) -> tuple[int, int]:
    return imzml_x - 1, imzml_y - 1


def find_spectrum_index_by_imzml_coord(
    exp: pyopenms.MSExperiment, x: int, y: int, z: int = 1
) -> int:
    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if not spec.metaValueExists("imzml:x") or not spec.metaValueExists("imzml:y"):
            continue
        spec_z = int(spec.getMetaValue("imzml:z")) if spec.metaValueExists("imzml:z") else 1
        if int(spec.getMetaValue("imzml:x")) == x and int(spec.getMetaValue("imzml:y")) == y and spec_z == z:
            return i
    raise RuntimeError(f"No spectrum at imzML pixel ({x}, {y}, {z})")


def demo_mode1_full_load(path: str) -> pyopenms.MSExperiment:
    print("=== Mode 1: ImzMLFile.load -> MSExperiment ===")
    exp = pyopenms.MSExperiment()
    pyopenms.ImzMLFile().load(path, exp)
    print(f"spectra: {exp.getNrSpectra()}")
    if exp.metaValueExists("imzml:imaging_mode"):
        print(f"imaging_mode: {exp.getMetaValue('imzml:imaging_mode')}")
    if exp.metaValueExists("imzml:max_count_x"):
        print(
            f"grid: {exp.getMetaValue('imzml:max_count_x')} x "
            f"{exp.getMetaValue('imzml:max_count_y')}"
        )
    if exp.getNrSpectra() > 0:
        print("spectrum[0] peaks:")
        print_peaks(exp.getSpectrum(0))
    return exp


def demo_mode2_file_handler(path: str) -> None:
    print("\n=== Mode 2: FileHandler.loadExperiment ===")
    print(f"file type: {pyopenms.FileHandler.getType(path)}")
    exp = pyopenms.MSExperiment()
    pyopenms.FileHandler().loadExperiment(path, exp)
    print(f"spectra via FileHandler: {exp.getNrSpectra()}")


def demo_mode3_streaming(path: str) -> int:
    print("\n=== Mode 3: ImzMLFile.load -> IMSDataConsumer ===")

    class Consumer:
        def __init__(self):
            self.count = 0
            self.first_size = 0

        def setExpectedSize(self, num_specs, num_chromo):
            pass

        def setExperimentalSettings(self, exp):
            pass

        def consumeChromatogram(self, chromo):
            pass

        def consumeSpectrum(self, spec):
            if self.count == 0:
                self.first_size = spec.size()
            self.count += 1

    consumer = Consumer()
    pyopenms.ImzMLFile().load(path, consumer)
    print(f"spectra streamed: {consumer.count}")
    print(f"first spectrum peaks: {consumer.first_size}")
    return consumer.count


def demo_mode4_on_disc_index(path: str) -> pyopenms.OnDiscImzMLExperiment:
    print("\n=== Mode 4: OnDiscImzMLExperiment.open (index + .ibd) ===")
    od = pyopenms.OnDiscImzMLExperiment()
    od.open(path)
    print(f"grid: {od.gridWidth()} x {od.gridHeight()}, spectra: {od.getNrSpectra()}")
    return od


def demo_mode5_on_disc_random_access(
    od: pyopenms.OnDiscImzMLExperiment,
    px0: int,
    py0: int,
    pz0: int,
    px1: int,
    py1: int,
    pz1: int,
) -> None:
    print("\n=== Mode 5: OnDisc random access (imzML 1-based pixels) ===")
    if od.getNrSpectra() == 0:
        print("  (no spectra)")
        return

    s0 = od.getSpectrum(0)
    print(f"getSpectrum(0): {s0.size()} peaks")

    by_coord = od.getSpectrumAtCoord(px0, py0, pz0)
    print(f"pixel ({px0},{py0},{pz0}) peaks:")
    print_peaks(by_coord)

    s_xy = od.getSpectrumAtCoord(px1, py1, pz1)
    print(f"pixel ({px1},{py1},{pz1}) peaks: {s_xy.size()}")


def demo_mode6_ram_lookup(
    path: str,
    px0: int,
    py0: int,
    px1: int,
    py1: int,
) -> pyopenms.MSImagingExperiment:
    print("\n=== Mode 6: MSImagingExperiment (RAM lookup, 0-based geometry) ===")
    imaging = pyopenms.MSImagingExperiment()
    pyopenms.ImzMLFile().load(path, imaging)
    geom = imaging.getGeometry()
    print(
        f"geometry: {geom.getWidth()} x {geom.getHeight()}, "
        f"pixels: {imaging.getNumberOfPixels()}"
    )
    gx0, gy0 = _geometry_coords(px0, py0)
    gx1, gy1 = _geometry_coords(px1, py1)
    if imaging.hasPixel(gx0, gy0):
        print(f"geometry pixel ({gx0},{gy0}) == imzML ({px0},{py0}) peaks:")
        print_peaks(imaging.getSpectrum(gx0, gy0))
    if imaging.hasPixel(gx1, gy1):
        print(f"geometry pixel ({gx1},{gy1}) == imzML ({px1},{py1}) peaks:")
        print_peaks(imaging.getSpectrum(gx1, gy1))
    return imaging


def demo_cross_mode_consistency(
    exp: pyopenms.MSExperiment,
    od: pyopenms.OnDiscImzMLExperiment,
    imaging: pyopenms.MSImagingExperiment,
    px: int,
    py: int,
    pz: int,
) -> None:
    print(f"\n=== Cross-mode check: pixel imzML ({px},{py},{pz}) ===")
    idx = find_spectrum_index_by_imzml_coord(exp, px, py, pz)
    full_size = exp.getSpectrum(idx).size()
    gx, gy = _geometry_coords(px, py)
    ram_size = imaging.getSpectrum(gx, gy).size()
    disc_size = od.getSpectrumAtCoord(px, py, pz).size()
    print(f"  MSExperiment index {idx}: {full_size} peaks")
    print(f"  MSImagingExperiment ({gx},{gy}): {ram_size} peaks")
    print(f"  OnDiscImzMLExperiment ({px},{py},{pz}): {disc_size} peaks")
    if full_size == ram_size == disc_size:
        print("  OK: peak counts match across all three access paths")
    else:
        print("  MISMATCH: peak counts differ", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Demonstrate imzML load modes in pyOpenMS.",
    )
    parser.add_argument(
        "imzml",
        nargs="?",
        default=_default_continuous_path(),
        help="Path to .imzML file (default: Example_Continuous test data)",
    )
    args = parser.parse_args()

    path = os.path.abspath(args.imzml)
    if not os.path.isfile(path):
        print(f"File not found: {path}", file=sys.stderr)
        sys.exit(1)

    ibd_path = _infer_ibd_path(path)
    if not os.path.isfile(ibd_path):
        print(f"Companion .ibd not found: {ibd_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading: {path}\n")

    exp = demo_mode1_full_load(path)
    demo_mode2_file_handler(path)
    demo_mode3_streaming(path)
    od = demo_mode4_on_disc_index(path)

    if od.getNrSpectra() == 0:
        print("No spectra in file; skipping coordinate-based demos.", file=sys.stderr)
        sys.exit(1)

    px0, py0, pz0 = _imzml_coords_from_index(od, 0)
    second_index = 1 if od.getNrSpectra() > 1 else 0
    px1, py1, pz1 = _imzml_coords_from_index(od, second_index)

    demo_mode5_on_disc_random_access(od, px0, py0, pz0, px1, py1, pz1)
    imaging = demo_mode6_ram_lookup(path, px0, py0, px1, py1)
    demo_cross_mode_consistency(exp, od, imaging, px1, py1, pz1)

    print("\nDone.")


if __name__ == "__main__":
    main()
