#!/usr/bin/env python3
"""
Comprehensive pyOpenMS benchmark suite.

Compares performance-critical paths across the Python bindings.
Benchmarks the nanobind-based Python bindings.

Usage:
    PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py
    PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --json results.json
    PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --quick
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import statistics
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Benchmark infrastructure
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkResult:
    name: str
    category: str
    iterations: int
    times_s: List[float] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def median_s(self) -> float:
        return statistics.median(self.times_s) if self.times_s else 0.0

    @property
    def mean_s(self) -> float:
        return statistics.mean(self.times_s) if self.times_s else 0.0

    @property
    def stdev_s(self) -> float:
        return statistics.stdev(self.times_s) if len(self.times_s) > 1 else 0.0

    @property
    def min_s(self) -> float:
        return min(self.times_s) if self.times_s else 0.0

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "category": self.category,
            "iterations": self.iterations,
            "median_s": self.median_s,
            "mean_s": self.mean_s,
            "stdev_s": self.stdev_s,
            "min_s": self.min_s,
            "times_s": self.times_s,
            "extra": self.extra,
        }


class BenchmarkSuite:
    def __init__(self, warmup: int = 1, repeats: int = 5, quick: bool = False):
        self.warmup = warmup
        self.repeats = repeats if not quick else 2
        self.quick = quick
        self.results: List[BenchmarkResult] = []

    def bench(self, name: str, category: str, fn: Callable, iterations: int = 1,
              setup: Optional[Callable] = None, teardown: Optional[Callable] = None,
              extra: Optional[Dict[str, Any]] = None) -> BenchmarkResult:
        """Run a benchmark function `repeats` times after `warmup` warmups."""
        result = BenchmarkResult(name=name, category=category, iterations=iterations,
                                 extra=extra or {})

        # Warmup
        for _ in range(self.warmup):
            if setup:
                setup()
            fn()
            if teardown:
                teardown()

        # Timed runs
        for _ in range(self.repeats):
            if setup:
                setup()
            gc.disable()
            t0 = time.perf_counter()
            for _ in range(iterations):
                fn()
            t1 = time.perf_counter()
            gc.enable()
            if teardown:
                teardown()
            result.times_s.append((t1 - t0) / iterations)

        self.results.append(result)
        med = result.median_s
        unit, scale = ("ms", 1e3) if med < 1.0 else ("s", 1.0)
        print(f"  {name:.<60s} {med * scale:>10.3f} {unit}  (stdev {result.stdev_s * scale:.3f} {unit})")
        return result

    def print_summary(self):
        print("\n" + "=" * 80)
        print("BENCHMARK SUMMARY")
        print("=" * 80)
        categories = sorted(set(r.category for r in self.results))
        for cat in categories:
            print(f"\n--- {cat} ---")
            for r in self.results:
                if r.category != cat:
                    continue
                med = r.median_s
                unit, scale = ("ms", 1e3) if med < 1.0 else ("s", 1.0)
                print(f"  {r.name:.<60s} {med * scale:>10.3f} {unit}")

    def save_json(self, path: str):
        data = {
            "metadata": {
                "python_version": sys.version,
                "numpy_version": np.__version__,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "repeats": self.repeats,
                "warmup": self.warmup,
            },
            "results": [r.to_dict() for r in self.results],
        }
        try:
            import pyopenms
            data["metadata"]["pyopenms_binding"] = getattr(pyopenms, "__binding__", "unknown")
        except Exception:
            pass
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        print(f"\nResults saved to {path}")


# ---------------------------------------------------------------------------
# Test data helpers
# ---------------------------------------------------------------------------

def find_mzml_file_auto() -> Optional[str]:
    """Find a suitable mzML file for benchmarking by searching known locations."""
    script_dir = Path(__file__).resolve().parent
    # Root of OpenMS repo (go up from src/pyOpenMS/tests/)
    repo_root = script_dir.parents[2]

    # Search order: prefer larger files for meaningful I/O benchmarks
    search_paths = [
        # pyOpenMS test files (always present)
        script_dir / "unittests" / "test2.mzML",           # 422 KB
        script_dir / "unittests" / "test.indexed.mzML",     # 640 KB
        # Larger TOPP test files
        repo_root / "src" / "tests" / "topp" / "FeatureFinderMetabo_1_input.mzML",
        repo_root / "src" / "tests" / "topp" / "OpenPepXL_input.mzML",
        # Example BSA data
        repo_root / "share" / "OpenMS" / "examples" / "BSA" / "BSA1.mzML",
        repo_root / "share" / "OpenMS" / "examples" / "ID" / "Ecoli_MS2_small.mzML",
    ]

    for p in search_paths:
        if p.exists() and p.stat().st_size > 1000:
            return str(p)

    # Fall back: search for any mzML > 10 KB
    for search_dir in [script_dir / "unittests", repo_root / "src" / "tests" / "topp"]:
        if search_dir.is_dir():
            for f in sorted(search_dir.glob("*.mzML"),
                            key=lambda p: p.stat().st_size, reverse=True):
                if f.stat().st_size > 10_000:
                    return str(f)
    return None


def build_synthetic_experiment(n_spectra: int = 500, peaks_per_spectrum: int = 1000,
                                n_chromatograms: int = 50, points_per_chrom: int = 200,
                                include_ms2: bool = True, include_im: bool = False):
    """Build a synthetic MSExperiment for benchmarking."""
    import pyopenms

    exp = pyopenms.MSExperiment()
    rng = np.random.default_rng(42)

    n_ms1 = n_spectra // 2 if include_ms2 else n_spectra
    n_ms2 = n_spectra - n_ms1 if include_ms2 else 0

    for i in range(n_ms1):
        spec = pyopenms.MSSpectrum()
        spec.setRT(float(i) * 0.5)
        spec.setMSLevel(1)
        spec.setNativeID(f"scan={i}")
        mz = np.sort(rng.uniform(100, 2000, peaks_per_spectrum))
        intensity = rng.exponential(1000, peaks_per_spectrum).astype(np.float32)
        spec.set_peaks((mz, intensity))

        if include_im:
            fda = pyopenms.FloatDataArray()
            fda.setName("Ion Mobility")
            im_values = rng.uniform(0.5, 1.5, peaks_per_spectrum).astype(np.float32)
            fda.set_data(im_values)
            spec.setFloatDataArrays([fda])

        exp.addSpectrum(spec)

    for i in range(n_ms2):
        spec = pyopenms.MSSpectrum()
        spec.setRT(float(i) * 0.5 + 0.25)
        spec.setMSLevel(2)
        spec.setNativeID(f"scan={n_ms1 + i}")
        mz = np.sort(rng.uniform(50, 1500, peaks_per_spectrum // 2))
        intensity = rng.exponential(500, peaks_per_spectrum // 2).astype(np.float32)
        spec.set_peaks((mz, intensity))

        prec = pyopenms.Precursor()
        prec.setMZ(rng.uniform(400, 1200))
        prec.setCharge(int(rng.choice([2, 3, 4])))
        prec.setIntensity(float(rng.exponential(10000)))
        spec.setPrecursors([prec])
        exp.addSpectrum(spec)

    for i in range(n_chromatograms):
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID(f"chrom_{i}")
        rts = np.linspace(0, 300, points_per_chrom)
        intensities = rng.exponential(100, points_per_chrom).astype(np.float64)
        chrom.set_peaks((rts, intensities))

        prec = pyopenms.Precursor()
        prec.setMZ(400 + i * 5.0)
        chrom.setPrecursor(prec)
        prod = pyopenms.Product()
        prod.setMZ(200 + i * 2.0)
        chrom.setProduct(prod)
        exp.addChromatogram(chrom)

    exp.updateRanges()
    return exp


def build_large_spectrum(n_peaks: int = 100_000):
    """Build a single large MSSpectrum."""
    import pyopenms

    rng = np.random.default_rng(123)
    spec = pyopenms.MSSpectrum()
    spec.setRT(42.0)
    spec.setMSLevel(1)
    spec.setNativeID("scan=1")
    mz = np.sort(rng.uniform(100, 2000, n_peaks))
    intensity = rng.exponential(1000, n_peaks).astype(np.float32)
    spec.set_peaks((mz, intensity))
    return spec


# ---------------------------------------------------------------------------
# Benchmark categories
# ---------------------------------------------------------------------------

def bench_file_io(suite: BenchmarkSuite, mzml_path: Optional[str]):
    """Benchmark file I/O operations."""
    import pyopenms
    print("\n[FILE I/O]")

    if mzml_path is None:
        print("  (skipped - no mzML test file found)")
        return

    file_size_mb = os.path.getsize(mzml_path) / (1024 * 1024)

    # MzMLFile load
    def load_mzml():
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(mzml_path, exp)
        return exp

    suite.bench("MzMLFile.load()", "File I/O", load_mzml,
                extra={"file": os.path.basename(mzml_path), "size_mb": round(file_size_mb, 1)})

    # Load once for write benchmark
    exp = load_mzml()

    # MzMLFile store
    tmp = tempfile.NamedTemporaryFile(suffix=".mzML", delete=False)
    tmp.close()
    try:
        suite.bench("MzMLFile.store()", "File I/O",
                    lambda: pyopenms.MzMLFile().store(tmp.name, exp),
                    extra={"n_spectra": exp.getNrSpectra()})
    finally:
        os.unlink(tmp.name)

    # IndexedMzMLFileLoader if indexed file available
    script_dir = Path(__file__).resolve().parent
    indexed_path = None
    candidate = script_dir / "unittests" / "test.indexed.mzML"
    if candidate.exists():
        indexed_path = str(candidate)

    if indexed_path:
        def load_indexed():
            exp2 = pyopenms.OnDiscMSExperiment()
            pyopenms.IndexedMzMLFileLoader().load(indexed_path, exp2)

        suite.bench("IndexedMzMLFileLoader.load()", "File I/O", load_indexed,
                    extra={"file": os.path.basename(indexed_path)})


def bench_spectrum_access(suite: BenchmarkSuite, exp, large_spec):
    """Benchmark spectrum data access patterns."""
    import pyopenms
    print("\n[SPECTRUM ACCESS]")

    n_spectra = exp.getNrSpectra()

    # --- get_peaks: the primary zero-copy path ---
    spec0 = exp[0]

    suite.bench(f"MSSpectrum.get_peaks() [{len(spec0)} peaks]", "Spectrum Access",
                lambda: spec0.get_peaks(), iterations=10_000)

    suite.bench(f"MSSpectrum.get_peaks() [{len(large_spec)} peaks]", "Spectrum Access",
                lambda: large_spec.get_peaks(), iterations=1000)

    # --- get_mz_array / get_intensity_array (individual) ---
    suite.bench(f"MSSpectrum.get_mz_array() [{len(large_spec)} peaks]", "Spectrum Access",
                lambda: large_spec.get_mz_array(), iterations=1000)

    suite.bench(f"MSSpectrum.get_intensity_array() [{len(large_spec)} peaks]", "Spectrum Access",
                lambda: large_spec.get_intensity_array(), iterations=1000)

    # --- Verify zero-copy: numpy operations on returned arrays ---
    def get_peaks_and_sum():
        mz, intensity = large_spec.get_peaks()
        return np.sum(mz) + np.sum(intensity)

    suite.bench(f"get_peaks() + np.sum [{len(large_spec)} peaks]", "Spectrum Access",
                get_peaks_and_sum, iterations=1000)

    # --- Verify array ownership (capsule): store and reuse ---
    def get_peaks_and_slice():
        mz, intensity = large_spec.get_peaks()
        return mz[100:200], intensity[100:200]

    suite.bench(f"get_peaks() + slice [{len(large_spec)} peaks]", "Spectrum Access",
                get_peaks_and_slice, iterations=1000)

    # --- Iteration over spectra (Python loop) ---
    def iterate_spectra_get_peaks():
        for spec in exp:
            mz, inten = spec.get_peaks()

    suite.bench(f"iterate all spectra + get_peaks() [n={n_spectra}]", "Spectrum Access",
                iterate_spectra_get_peaks)

    # --- __getitem__(index) vs getSpectra()[index] ---
    suite.bench(f"MSExperiment[i] x {n_spectra}", "Spectrum Access",
                lambda: [exp[i] for i in range(n_spectra)])

    def get_spectra_vector_index():
        spectra = exp.getSpectra()
        for i in range(len(spectra)):
            _ = spectra[i]

    suite.bench(f"getSpectra() + index [n={n_spectra}]", "Spectrum Access",
                get_spectra_vector_index)

    # --- __getitem__ on spectrum (peak-by-peak, slow path) ---
    def getitem_loop():
        for i in range(min(10_000, len(large_spec))):
            _ = large_spec[i]

    suite.bench("MSSpectrum.__getitem__ x 10000 peaks", "Spectrum Access",
                getitem_loop, iterations=10)

    # --- Metadata access ---
    def metadata_access():
        for spec in exp:
            _ = spec.getRT()
            _ = spec.getMSLevel()
            _ = spec.getNativeID()

    suite.bench(f"RT/MSLevel/NativeID access [n={n_spectra}]", "Spectrum Access",
                metadata_access, iterations=100)


def bench_spectrum_mutation(suite: BenchmarkSuite):
    """Benchmark spectrum creation and mutation."""
    import pyopenms
    print("\n[SPECTRUM MUTATION]")

    rng = np.random.default_rng(99)
    n_peaks = 10_000

    # --- set_peaks from tuple of numpy arrays ---
    mz = np.sort(rng.uniform(100, 2000, n_peaks))
    intensity = rng.exponential(1000, n_peaks).astype(np.float32)

    def set_peaks_tuple():
        spec = pyopenms.MSSpectrum()
        spec.set_peaks((mz, intensity))

    suite.bench(f"set_peaks((mz, int)) [n={n_peaks}]", "Spectrum Mutation",
                set_peaks_tuple, iterations=100)

    # --- push_back (peak by peak) ---
    def push_back_loop():
        spec = pyopenms.MSSpectrum()
        for i in range(10_000):
            p = pyopenms.Peak1D()
            p.setMZ(float(i))
            p.setIntensity(float(i))
            spec.push_back(p)

    suite.bench("push_back(Peak1D) x 10000", "Spectrum Mutation",
                push_back_loop, iterations=10)

    # --- MSSpectrum creation ---
    suite.bench("MSSpectrum() construction x 10000", "Spectrum Mutation",
                lambda: [pyopenms.MSSpectrum() for _ in range(10_000)], iterations=10)

    # --- Peak1D creation ---
    suite.bench("Peak1D() construction x 10000", "Spectrum Mutation",
                lambda: [pyopenms.Peak1D() for _ in range(10_000)], iterations=10)

    # --- sortByPosition ---
    spec_unsorted = pyopenms.MSSpectrum()
    mz_unsorted = rng.uniform(100, 2000, n_peaks)
    int_unsorted = rng.exponential(1000, n_peaks).astype(np.float32)

    def sort_spectrum():
        spec_unsorted.set_peaks((mz_unsorted, int_unsorted))
        spec_unsorted.sortByPosition()

    suite.bench(f"set_peaks + sortByPosition [n={n_peaks}]", "Spectrum Mutation",
                sort_spectrum, iterations=50)

    # --- findNearest ---
    spec_sorted = pyopenms.MSSpectrum()
    spec_sorted.set_peaks((np.sort(mz_unsorted), int_unsorted))
    targets = rng.uniform(200, 1800, 1000)

    def find_nearest_loop():
        for t in targets:
            spec_sorted.findNearest(float(t))

    suite.bench("findNearest x 1000", "Spectrum Mutation",
                find_nearest_loop, iterations=10)

    # --- calculateTIC ---
    large = build_large_spectrum(50_000)
    suite.bench("calculateTIC [50k peaks]", "Spectrum Mutation",
                lambda: large.calculateTIC(), iterations=100)


def bench_chromatogram_access(suite: BenchmarkSuite, exp):
    """Benchmark chromatogram access patterns."""
    import pyopenms
    print("\n[CHROMATOGRAM ACCESS]")

    n_chrom = exp.getNrChromatograms()
    if n_chrom == 0:
        print("  (skipped - no chromatograms)")
        return

    # --- get_peaks ---
    chrom0 = exp.getChromatogram(0)

    suite.bench(f"MSChromatogram.get_peaks() [{len(chrom0)} pts]", "Chromatogram Access",
                lambda: chrom0.get_peaks(), iterations=10_000)

    # --- Iterate all chromatograms ---
    def iterate_chroms():
        for i in range(n_chrom):
            chrom = exp.getChromatogram(i)
            rt, inten = chrom.get_peaks()

    suite.bench(f"iterate chromatograms + get_peaks [n={n_chrom}]", "Chromatogram Access",
                iterate_chroms, iterations=100)

    # --- Chromatogram set_peaks ---
    rng = np.random.default_rng(77)
    rt_arr = np.linspace(0, 600, 5000)
    int_arr = rng.exponential(100, 5000)

    def chrom_set_peaks():
        c = pyopenms.MSChromatogram()
        c.set_peaks((rt_arr, int_arr))

    suite.bench("MSChromatogram.set_peaks [5000 pts]", "Chromatogram Access",
                chrom_set_peaks, iterations=100)


def bench_data_arrays(suite: BenchmarkSuite):
    """Benchmark DataArray operations (FloatDataArray, IntegerDataArray, StringDataArray)."""
    import pyopenms
    print("\n[DATA ARRAYS]")

    n = 100_000
    rng = np.random.default_rng(55)

    # --- FloatDataArray ---
    fda = pyopenms.FloatDataArray()
    fda_init_data = rng.uniform(0, 100, n).astype(np.float32)
    fda.set_data(fda_init_data)

    suite.bench(f"FloatDataArray.get_data() [n={n}]", "Data Arrays",
                lambda: fda.get_data(), iterations=100)

    # --- FloatDataArray get_data_view (zero-copy memory view) ---
    if hasattr(fda, 'get_data_view'):
        suite.bench(f"FloatDataArray.get_data_view() [n={n}]", "Data Arrays",
                    lambda: fda.get_data_view(), iterations=10_000)

        # Compare: get_data (copy) vs get_data_view (zero-copy)
        def use_get_data():
            arr = np.array(fda.get_data(), dtype=np.float32)
            return np.sum(arr)

        def use_get_data_view():
            arr = fda.get_data_view()
            return np.sum(arr)

        suite.bench(f"FloatDataArray get_data+np.sum [n={n}]", "Data Arrays",
                    use_get_data, iterations=100)

        suite.bench(f"FloatDataArray get_data_view+np.sum [n={n}]", "Data Arrays",
                    use_get_data_view, iterations=100)

    # --- FloatDataArray set_data ---
    data_np = rng.uniform(0, 100, n).astype(np.float32)
    suite.bench(f"FloatDataArray.set_data() [n={n}]", "Data Arrays",
                lambda: fda.set_data(data_np), iterations=100)

    # --- FloatDataArray element-by-element access ---
    def fda_getitem_loop():
        for i in range(min(10_000, len(fda))):
            _ = fda[i]

    suite.bench("FloatDataArray.__getitem__ x 10000", "Data Arrays",
                fda_getitem_loop, iterations=10)

    def fda_setitem_loop():
        for i in range(min(10_000, len(fda))):
            fda[i] = float(i)

    suite.bench("FloatDataArray.__setitem__ x 10000", "Data Arrays",
                fda_setitem_loop, iterations=10)

    # --- IntegerDataArray ---
    ida = pyopenms.IntegerDataArray()
    ida_data = np.arange(n, dtype=np.int32)
    ida.set_data(ida_data)

    suite.bench(f"IntegerDataArray.get_data() [n={n}]", "Data Arrays",
                lambda: ida.get_data(), iterations=100)

    suite.bench(f"IntegerDataArray.set_data() [n={n}]", "Data Arrays",
                lambda: ida.set_data(ida_data), iterations=100)

    # --- StringDataArray ---
    sda = pyopenms.StringDataArray()
    n_strings = min(n, 10_000)
    strings = [f"ion_{i}" for i in range(n_strings)]
    if hasattr(sda, 'set_data'):
        sda.set_data(strings)
    else:
        for s in strings:
            sda.push_back(s)

    suite.bench(f"StringDataArray.get_data() [n={n_strings}]", "Data Arrays",
                lambda: sda.get_data(), iterations=100)

    if hasattr(sda, 'set_data'):
        suite.bench(f"StringDataArray.set_data() [n={n_strings}]", "Data Arrays",
                    lambda: sda.set_data(strings), iterations=100)


def bench_memory_views(suite: BenchmarkSuite):
    """Benchmark zero-copy memory views vs copy-based access."""
    import pyopenms
    print("\n[MEMORY VIEWS / ZERO-COPY]")

    n_peaks = 100_000
    rng = np.random.default_rng(42)

    # Build spectrum with ion mobility data
    spec = pyopenms.MSSpectrum()
    spec.setRT(100.0)
    spec.setMSLevel(1)
    mz = np.sort(rng.uniform(100, 2000, n_peaks))
    intensity = rng.exponential(1000, n_peaks).astype(np.float32)
    spec.set_peaks((mz, intensity))

    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    im_data = rng.uniform(0.5, 1.5, n_peaks).astype(np.float32)
    fda.set_data(im_data)
    spec.setFloatDataArrays([fda])

    # --- get_peaks (capsule-based ownership) ---
    suite.bench(f"get_peaks() zero-copy [n={n_peaks}]", "Memory Views",
                lambda: spec.get_peaks(), iterations=1000)

    # --- get_mz_array vs get_peaks (compare overhead) ---
    suite.bench(f"get_mz_array() [n={n_peaks}]", "Memory Views",
                lambda: spec.get_mz_array(), iterations=1000)

    suite.bench(f"get_intensity_array() [n={n_peaks}]", "Memory Views",
                lambda: spec.get_intensity_array(), iterations=1000)

    # --- get_drift_time_array (copy) vs get_drift_time_array_view (zero-copy view) ---
    if hasattr(spec, 'get_drift_time_array'):
        suite.bench(f"get_drift_time_array() copy [n={n_peaks}]", "Memory Views",
                    lambda: spec.get_drift_time_array(), iterations=1000)

    if hasattr(spec, 'get_drift_time_array_view'):
        suite.bench(f"get_drift_time_array_view() zero-copy [n={n_peaks}]", "Memory Views",
                    lambda: spec.get_drift_time_array_view(), iterations=10_000)

    # --- FloatDataArray: get_data (copy) vs get_data_view (zero-copy view) ---
    fda_large = pyopenms.FloatDataArray()
    large_data = rng.uniform(0, 100, n_peaks).astype(np.float32)
    fda_large.set_data(large_data)

    suite.bench(f"FDA.get_data() copy [n={n_peaks}]", "Memory Views",
                lambda: fda_large.get_data(), iterations=100)

    if hasattr(fda_large, 'get_data_view'):
        suite.bench(f"FDA.get_data_view() zero-copy [n={n_peaks}]", "Memory Views",
                    lambda: fda_large.get_data_view(), iterations=10_000)

    # --- Numpy operations on zero-copy vs copy ---
    def copy_then_compute():
        arr = np.array(fda_large.get_data(), dtype=np.float32)
        return float(np.mean(arr))

    def zerocopy_then_compute():
        arr = fda_large.get_data_view()
        return float(np.mean(arr))

    suite.bench(f"copy + np.mean [n={n_peaks}]", "Memory Views",
                copy_then_compute, iterations=50)

    if hasattr(fda_large, 'get_data_view'):
        suite.bench(f"zero-copy + np.mean [n={n_peaks}]", "Memory Views",
                    zerocopy_then_compute, iterations=50)

    # --- Verify numpy array flags ---
    mz_arr, int_arr = spec.get_peaks()
    info = {
        "mz_contiguous": bool(mz_arr.flags['C_CONTIGUOUS']),
        "mz_writable": bool(mz_arr.flags['WRITEABLE']),
        "mz_owndata": bool(mz_arr.flags['OWNDATA']),
        "mz_dtype": str(mz_arr.dtype),
        "int_dtype": str(int_arr.dtype),
    }
    print(f"  [info] get_peaks array flags: {info}")

    if hasattr(fda_large, 'get_data_view'):
        mv = fda_large.get_data_view()
        mv_info = {
            "contiguous": bool(mv.flags['C_CONTIGUOUS']),
            "writable": bool(mv.flags['WRITEABLE']),
            "dtype": str(mv.dtype),
        }
        print(f"  [info] get_data_view array flags: {mv_info}")

    # --- Round-trip: get_peaks -> modify numpy -> set_peaks ---
    def roundtrip():
        mz, inten = spec.get_peaks()
        mz = mz * 1.001  # creates a new array
        inten = inten * 1.1
        spec.set_peaks((mz, inten))

    suite.bench(f"get_peaks + modify + set_peaks [n={n_peaks}]", "Memory Views",
                roundtrip, iterations=50)


def bench_arrow_export(suite: BenchmarkSuite, exp):
    """Benchmark Arrow zero-copy export."""
    print("\n[ARROW EXPORT]")

    try:
        import pyarrow as pa
    except ImportError:
        print("  (skipped - pyarrow not installed)")
        return

    n_spectra = exp.getNrSpectra()
    n_chroms = exp.getNrChromatograms()

    # --- Check if zero-copy path is available AND functional ---
    has_zerocopy = False
    try:
        from pyopenms._arrow_zerocopy import spectra_to_arrow, chromatograms_to_arrow
        # Verify it actually works
        spectra_to_arrow(exp, format='long', ms_levels=[1])
        has_zerocopy = True
        print("  [info] Zero-copy C++ Arrow export available and functional")
    except (ImportError, RuntimeError):
        has_zerocopy = False
        print("  [info] Using Python fallback Arrow export")
        # If _arrow_zerocopy import succeeds but cast fails, block it so to_arrow
        # uses the Python fallback path
        import sys
        sys.modules['pyopenms._arrow_zerocopy'] = None

    # --- Spectra to Arrow: long format ---
    suite.bench(f"to_arrow(spectra, long) [n={n_spectra}]", "Arrow Export",
                lambda: exp.to_arrow(data='spectra', format='long'))

    # --- Spectra to Arrow: semi-wide format ---
    suite.bench(f"to_arrow(spectra, semi_wide) [n={n_spectra}]", "Arrow Export",
                lambda: exp.to_arrow(data='spectra', format='semi_wide'))

    # --- With column filtering ---
    suite.bench(f"to_arrow(spectra, long, cols=[mz,int,rt]) [n={n_spectra}]", "Arrow Export",
                lambda: exp.to_arrow(data='spectra', format='long',
                                     columns=['mz', 'intensity', 'rt']))

    # --- With MS level filtering ---
    suite.bench(f"to_arrow(spectra, long, ms_levels=[1]) [n={n_spectra}]", "Arrow Export",
                lambda: exp.to_arrow(data='spectra', format='long', ms_levels=[1]))

    # --- With RT range filtering ---
    exp.updateRanges()
    mid_rt = (exp.getMinRT() + exp.getMaxRT()) / 2
    suite.bench(f"to_arrow(spectra, long, rt_range) [n={n_spectra}]", "Arrow Export",
                lambda: exp.to_arrow(data='spectra', format='long',
                                     min_rt=exp.getMinRT(), max_rt=mid_rt))

    # --- Without precursor info (less work) ---
    suite.bench(f"to_arrow(spectra, long, no_prec) [n={n_spectra}]", "Arrow Export",
                lambda: exp.to_arrow(data='spectra', format='long',
                                     include_precursor_info=False))

    # --- Chromatograms to Arrow ---
    if n_chroms > 0:
        suite.bench(f"to_arrow(chroms, long) [n={n_chroms}]", "Arrow Export",
                    lambda: exp.to_arrow(data='chromatograms', format='long'))

        suite.bench(f"to_arrow(chroms, semi_wide) [n={n_chroms}]", "Arrow Export",
                    lambda: exp.to_arrow(data='chromatograms', format='semi_wide'))

    # --- Both spectra and chromatograms ---
    suite.bench(f"to_arrow(both, long)", "Arrow Export",
                lambda: exp.to_arrow(data='both', format='long'))

    # --- If zero-copy available, also bench the raw C++ function ---
    if has_zerocopy:
        suite.bench(f"spectra_to_arrow() C++ direct [n={n_spectra}]", "Arrow Export",
                    lambda: spectra_to_arrow(exp, format='long'))

        if n_chroms > 0:
            suite.bench(f"chromatograms_to_arrow() C++ direct [n={n_chroms}]", "Arrow Export",
                        lambda: chromatograms_to_arrow(exp, format='long'))

    # --- Arrow table operations (measure downstream usability) ---
    tbl = exp.to_arrow(data='spectra', format='long')
    suite.bench("Arrow Table .to_pandas()", "Arrow Export",
                lambda: tbl.to_pandas(), iterations=3)

    suite.bench("Arrow Table .to_pydict()", "Arrow Export",
                lambda: tbl.to_pydict(), iterations=10)

    suite.bench("Arrow Table .column('mz').to_numpy()", "Arrow Export",
                lambda: tbl.column('mz').to_numpy(), iterations=100,
                extra={"n_rows": tbl.num_rows})

    # --- Single spectrum to Arrow ---
    import pyopenms
    spec = exp[0]
    suite.bench(f"MSSpectrum.to_arrow() [{len(spec)} peaks]", "Arrow Export",
                lambda: spec.to_arrow(), iterations=100)


def bench_dataframe_export(suite: BenchmarkSuite, exp, large_spec):
    """Benchmark DataFrame export."""
    print("\n[DATAFRAME EXPORT]")

    try:
        import pandas as pd
    except ImportError:
        print("  (skipped - pandas not installed)")
        return

    n_spectra = exp.getNrSpectra()

    # --- MSSpectrum.to_df ---
    suite.bench(f"MSSpectrum.to_df() [{len(large_spec)} peaks]", "DataFrame Export",
                lambda: large_spec.to_df(), iterations=20)

    suite.bench(f"MSSpectrum.to_df(columns=[mz,int]) [{len(large_spec)} peaks]", "DataFrame Export",
                lambda: large_spec.to_df(columns=['mz', 'intensity']), iterations=20)

    # --- MSExperiment.to_df (wide format: one row per spectrum) ---
    suite.bench(f"MSExperiment.to_df(wide) [n={n_spectra}]", "DataFrame Export",
                lambda: exp.to_df(long_format=False))

    # --- MSExperiment.to_df (long format: one row per peak) ---
    suite.bench(f"MSExperiment.to_df(long) [n={n_spectra}]", "DataFrame Export",
                lambda: exp.to_df(long_format=True))

    # --- MSExperiment.to_df with ms_level filter ---
    suite.bench(f"MSExperiment.to_df(long, ms1 only) [n={n_spectra}]", "DataFrame Export",
                lambda: exp.to_df(long_format=True, ms_levels=[1]))

    # --- MSChromatogram.to_df ---
    import pyopenms
    if exp.getNrChromatograms() > 0:
        chrom = exp.getChromatogram(0) if hasattr(exp, 'getChromatogram') else exp.getChromatograms()[0]
        suite.bench(f"MSChromatogram.to_df() [{len(chrom)} pts]", "DataFrame Export",
                    lambda: chrom.to_df(), iterations=100)

    # --- get2DPeakDataLong ---
    exp.updateRanges()
    suite.bench(f"get2DPeakDataLong (MS1) [n={n_spectra}]", "DataFrame Export",
                lambda: exp.get2DPeakDataLong(
                    exp.getMinRT(), exp.getMaxRT(),
                    exp.getMinMZ(), exp.getMaxMZ(), 1))

    # --- get_data_dict (intermediate step) ---
    suite.bench(f"MSSpectrum.get_data_dict() [{len(large_spec)} peaks]", "DataFrame Export",
                lambda: large_spec.get_data_dict(), iterations=20)


def bench_rasterize(suite: BenchmarkSuite, exp):
    """Benchmark 2D rasterization."""
    import pyopenms
    print("\n[RASTERIZE]")

    exp.updateRanges()
    if exp.getNrSpectra() == 0:
        print("  (skipped - no spectra)")
        return

    # Count total MS1 peaks
    total_peaks = sum(len(s) for s in exp if s.getMSLevel() == 1)

    min_rt, max_rt = exp.getMinRT(), exp.getMaxRT()
    min_mz, max_mz = exp.getMinMZ(), exp.getMaxMZ()

    print(f"  [info] {total_peaks} total MS1 peaks across {exp.getNrSpectra()} spectra")

    # --- Small grid ---
    output_small = np.zeros((100, 100), dtype=np.float32)
    suite.bench(f"rasterizeRTMZ 100x100 sum [{total_peaks} peaks]", "Rasterize",
                lambda: exp.rasterizeRTMZ(output_small, min_rt, max_rt,
                                           min_mz, max_mz, 1, "sum"),
                iterations=10)

    # --- Medium grid ---
    output_med = np.zeros((500, 500), dtype=np.float32)
    suite.bench(f"rasterizeRTMZ 500x500 sum [{total_peaks} peaks]", "Rasterize",
                lambda: exp.rasterizeRTMZ(output_med, min_rt, max_rt,
                                           min_mz, max_mz, 1, "sum"),
                iterations=5)

    # --- Large grid ---
    output_large = np.zeros((1000, 1000), dtype=np.float32)
    suite.bench(f"rasterizeRTMZ 1000x1000 sum [{total_peaks} peaks]", "Rasterize",
                lambda: exp.rasterizeRTMZ(output_large, min_rt, max_rt,
                                           min_mz, max_mz, 1, "sum"),
                iterations=3)

    # --- Extra large grid ---
    output_xlarge = np.zeros((2500, 2500), dtype=np.float32)
    suite.bench(f"rasterizeRTMZ 2500x2500 sum [{total_peaks} peaks]", "Rasterize",
                lambda: exp.rasterizeRTMZ(output_xlarge, min_rt, max_rt,
                                           min_mz, max_mz, 1, "sum"),
                iterations=3)

    # --- Max aggregation ---
    suite.bench(f"rasterizeRTMZ 500x500 max [{total_peaks} peaks]", "Rasterize",
                lambda: exp.rasterizeRTMZ(output_med, min_rt, max_rt,
                                           min_mz, max_mz, 1, "max"),
                iterations=5)


def bench_experiment_operations(suite: BenchmarkSuite, exp):
    """Benchmark MSExperiment-level operations."""
    import pyopenms
    print("\n[EXPERIMENT OPERATIONS]")

    n_spectra = exp.getNrSpectra()

    # --- Copy experiment ---
    suite.bench(f"MSExperiment copy [n={n_spectra}]", "Experiment Ops",
                lambda: pyopenms.MSExperiment(exp), iterations=3)

    # --- updateRanges ---
    suite.bench(f"MSExperiment.updateRanges() [n={n_spectra}]", "Experiment Ops",
                lambda: exp.updateRanges(), iterations=10)

    # --- sortSpectra ---
    def sort_copy():
        e = pyopenms.MSExperiment(exp)
        e.sortSpectra(True)

    suite.bench(f"sortSpectra(sort_mz=True) [n={n_spectra}]", "Experiment Ops",
                sort_copy, iterations=3)

    # --- getMSLevels ---
    suite.bench(f"getMSLevels() [n={n_spectra}]", "Experiment Ops",
                lambda: exp.getMSLevels(), iterations=10)

    # --- size / getNrSpectra ---
    suite.bench("getNrSpectra() x 100000", "Experiment Ops",
                lambda: [exp.getNrSpectra() for _ in range(100_000)])

    # --- Build experiment from scratch ---
    rng = np.random.default_rng(12)

    def build_experiment():
        e = pyopenms.MSExperiment()
        for i in range(100):
            s = pyopenms.MSSpectrum()
            s.setRT(float(i))
            s.setMSLevel(1)
            mz = np.sort(rng.uniform(100, 2000, 500))
            inten = rng.exponential(1000, 500).astype(np.float32)
            s.set_peaks((mz, inten))
            e.addSpectrum(s)
        return e

    suite.bench("Build MSExperiment (100 spectra x 500 peaks)", "Experiment Ops",
                build_experiment, iterations=5)


def bench_type_conversions(suite: BenchmarkSuite):
    """Benchmark type conversion overhead."""
    import pyopenms
    print("\n[TYPE CONVERSIONS]")

    # --- String conversion ---
    spec = pyopenms.MSSpectrum()
    spec.setNativeID("scan=12345")

    suite.bench("getNativeID() str conversion x 10000", "Type Conversions",
                lambda: [spec.getNativeID() for _ in range(10000)])

    suite.bench("setNativeID() str conversion x 10000", "Type Conversions",
                lambda: [spec.setNativeID(f"scan={i}") for i in range(10000)])

    # --- MetaValue get/set ---
    spec.setMetaValue("test_int", 42)
    spec.setMetaValue("test_float", 3.14)
    spec.setMetaValue("test_string", "hello world")

    suite.bench("getMetaValue(int) x 10000", "Type Conversions",
                lambda: [spec.getMetaValue("test_int") for _ in range(10000)])

    suite.bench("getMetaValue(float) x 10000", "Type Conversions",
                lambda: [spec.getMetaValue("test_float") for _ in range(10000)])

    suite.bench("getMetaValue(str) x 10000", "Type Conversions",
                lambda: [spec.getMetaValue("test_string") for _ in range(10000)])

    suite.bench("setMetaValue x 10000", "Type Conversions",
                lambda: [spec.setMetaValue("bench_key", i) for i in range(10000)])


def bench_precursor_access(suite: BenchmarkSuite, exp):
    """Benchmark precursor/product access patterns."""
    import pyopenms
    print("\n[PRECURSOR ACCESS]")

    # Find MS2 spectra
    ms2_specs = [s for s in exp if s.getMSLevel() == 2]
    if not ms2_specs:
        print("  (skipped - no MS2 spectra)")
        return

    n_ms2 = len(ms2_specs)

    # --- getPrecursors on single spectrum ---
    suite.bench(f"getPrecursors() x {n_ms2}", "Precursor Access",
                lambda: [s.getPrecursors() for s in ms2_specs], iterations=100)

    # --- Full precursor info extraction ---
    def extract_precursor_info():
        results = []
        for s in ms2_specs:
            precs = s.getPrecursors()
            if precs:
                p = precs[0]
                results.append((p.getMZ(), p.getCharge(), p.getIntensity()))
        return results

    suite.bench(f"extract precursor mz/charge/int x {n_ms2}", "Precursor Access",
                extract_precursor_info, iterations=100)

    # --- setPrecursors ---
    prec = pyopenms.Precursor()
    prec.setMZ(500.0)
    prec.setCharge(2)

    def set_precursors():
        for s in ms2_specs:
            s.setPrecursors([prec])

    suite.bench(f"setPrecursors x {n_ms2}", "Precursor Access",
                set_precursors, iterations=100)


def bench_identification(suite: BenchmarkSuite):
    """Benchmark PeptideIdentification/Hit operations."""
    import pyopenms
    print("\n[IDENTIFICATION]")

    # --- Build PeptideIdentificationList ---
    n_ids = 500

    def build_id_list():
        pil = pyopenms.PeptideIdentificationList()
        for i in range(n_ids):
            pid = pyopenms.PeptideIdentification()
            pid.setRT(float(i))
            pid.setMZ(400.0 + i * 0.5)
            hit = pyopenms.PeptideHit()
            hit.setScore(0.99 - i * 0.001)
            hit.setCharge(2)
            hit.setSequence(pyopenms.AASequence.fromString("PEPTIDER"))
            pid.setHits([hit])
            pil.append(pid)
        return pil

    suite.bench(f"Build PeptideIdentificationList [n={n_ids}]", "Identification",
                build_id_list, iterations=3)

    pil = build_id_list()

    # --- to_df ---
    if hasattr(pil, 'to_df'):
        suite.bench(f"PeptideIdentificationList.to_df() [n={n_ids}]", "Identification",
                    lambda: pil.to_df(), iterations=3)

    # --- AASequence operations ---
    suite.bench("AASequence.fromString x 1000", "Identification",
                lambda: [pyopenms.AASequence.fromString("PEPTIDER") for _ in range(1000)],
                iterations=5)

    seq = pyopenms.AASequence.fromString("PEPTIDER")
    suite.bench("AASequence.toString x 1000", "Identification",
                lambda: [seq.toString() for _ in range(1000)],
                iterations=5)

    suite.bench("AASequence.getMonoWeight x 10000", "Identification",
                lambda: [seq.getMonoWeight() for _ in range(10000)])


def bench_massql_df(suite: BenchmarkSuite, exp):
    """Benchmark MassQL DataFrame export."""
    print("\n[MASSQL DATAFRAME]")

    if not hasattr(exp, 'get_massql_df'):
        print("  (skipped - get_massql_df not available)")
        return

    n_spectra = exp.getNrSpectra()
    suite.bench(f"get_massql_df() [n={n_spectra}]", "MassQL DataFrame",
                lambda: exp.get_massql_df(), iterations=3)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="pyOpenMS Benchmark Suite")
    parser.add_argument("--json", type=str, help="Save results to JSON file")
    parser.add_argument("--quick", action="store_true", help="Reduced iterations for quick run")
    parser.add_argument("--repeats", type=int, default=5, help="Number of repeat measurements")
    parser.add_argument("--mzml", type=str, help="Path to mzML file for I/O benchmarks")
    parser.add_argument("--n-spectra", type=int, default=1000,
                        help="Number of spectra in synthetic experiment")
    parser.add_argument("--peaks-per-spectrum", type=int, default=2000,
                        help="Peaks per spectrum in synthetic experiment")
    parser.add_argument("--categories", type=str, nargs="+",
                        help="Only run specific categories (e.g., 'arrow' 'memory')")
    args = parser.parse_args()

    print("=" * 80)
    print("pyOpenMS Comprehensive Benchmark Suite")
    print("=" * 80)

    # Import and report binding info
    import pyopenms
    binding = getattr(pyopenms, "__binding__", "unknown")
    print(f"Binding type: {binding}")
    print(f"Python: {sys.version}")
    print(f"NumPy: {np.__version__}")
    try:
        import pyarrow
        print(f"PyArrow: {pyarrow.__version__}")
    except ImportError:
        print("PyArrow: not installed")
    try:
        import pandas
        print(f"Pandas: {pandas.__version__}")
    except ImportError:
        print("Pandas: not installed")

    # Check zero-copy availability
    try:
        from pyopenms._arrow_zerocopy import spectra_to_arrow
        print("Arrow zero-copy: available")
    except ImportError:
        print("Arrow zero-copy: NOT available")

    suite = BenchmarkSuite(warmup=1, repeats=args.repeats, quick=args.quick)

    # --- Determine mzML path ---
    mzml_path = args.mzml
    if not mzml_path:
        mzml_path = find_mzml_file_auto()
    if mzml_path:
        print(f"Test mzML: {mzml_path} ({os.path.getsize(mzml_path) / 1024:.0f} KB)")
    else:
        print("Test mzML: none found (I/O benchmarks will be skipped)")

    # --- Build synthetic data ---
    n_spec = args.n_spectra if not args.quick else 100
    pps = args.peaks_per_spectrum if not args.quick else 500
    print(f"\nBuilding synthetic experiment ({n_spec} spectra x {pps} peaks)...")
    exp = build_synthetic_experiment(
        n_spectra=n_spec,
        peaks_per_spectrum=pps,
        n_chromatograms=50 if not args.quick else 10,
        points_per_chrom=200,
        include_ms2=True,
        include_im=False,
    )
    print(f"  -> {exp.getNrSpectra()} spectra, {exp.getNrChromatograms()} chromatograms")

    large_n = 100_000 if not args.quick else 20_000
    print(f"Building large spectrum ({large_n} peaks)...")
    large_spec = build_large_spectrum(large_n)

    def should_run(cat: str) -> bool:
        if not args.categories:
            return True
        return any(c.lower() in cat.lower() for c in args.categories)

    # --- Run benchmarks ---
    print("\n" + "-" * 80)

    if should_run("file"):
        bench_file_io(suite, mzml_path)

    if should_run("spectrum"):
        bench_spectrum_access(suite, exp, large_spec)

    if should_run("mutation"):
        bench_spectrum_mutation(suite)

    if should_run("chromatogram"):
        bench_chromatogram_access(suite, exp)

    if should_run("data array"):
        bench_data_arrays(suite)

    if should_run("memory") or should_run("zero"):
        bench_memory_views(suite)

    if should_run("arrow"):
        bench_arrow_export(suite, exp)

    if should_run("dataframe") or should_run("df"):
        bench_dataframe_export(suite, exp, large_spec)

    if should_run("raster"):
        bench_rasterize(suite, exp)

    if should_run("experiment") or should_run("ops"):
        bench_experiment_operations(suite, exp)

    if should_run("type") or should_run("conversion"):
        bench_type_conversions(suite)

    if should_run("precursor"):
        bench_precursor_access(suite, exp)

    if should_run("ident"):
        bench_identification(suite)

    if should_run("massql"):
        bench_massql_df(suite, exp)

    # --- Summary ---
    suite.print_summary()

    if args.json:
        suite.save_json(args.json)

    return 0


if __name__ == "__main__":
    sys.exit(main())
