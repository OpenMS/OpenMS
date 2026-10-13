#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Benchmark the in-process vendor readers and prove that their output did not change.

The script serves three purposes for work on ``BrukerTimsFile`` (.d) and
``ThermoRawFile`` (.raw):

1. Benchmark: wall time, CPU time and peak memory of one reader call per
   reader option / export mode, each measured in a fresh process.
2. Output parity: an exact content fingerprint of what a reader returned
   (spectrum and peak counts, SHA-256 of the m/z, intensity, ion mobility and
   other data arrays, of the spectrum metadata and of the experiment settings),
   so that a report taken before a code change can be compared with one taken
   after it. ``compare`` says IDENTICAL, REORDERED (same content, different
   order of spectra or of peaks within spectra) or DIFFERENT, and names the
   component and the first spectra that differ.
3. Cross-reader comparison with tolerances (``diff``), for example the
   in-process Thermo reader against an mzML written by ThermoRawFileParser.

Typical workflow::

    BUILD=/path/to/OpenMS-build
    export PYTHONPATH=$BUILD/pyOpenMS      # or: pip install pyopenms
    DDA=/data/DDA_HeLa.d  DIA=/data/DIA_HeLa.d  RAW=/data/run.raw
    S=tools/scripts/benchmark_vendor_readers.py

    python3 $S list                                  # the case matrix

    # before touching the reader: the baseline
    python3 $S run --bruker-dda $DDA --bruker-dia $DIA --thermo $RAW \\
        --openms-bin $BUILD/bin --repeat 3 -o baseline.json

    # after each change: the quick parity loop (small frame windows + Thermo)
    python3 $S run --bruker-dda $DDA --bruker-dia $DIA --thermo $RAW \\
        --filter '@|^thermo/' --parity-only --jobs 4 -o quick.json
    python3 $S compare baseline.json quick.json --allow-missing

    # at the end: everything again, with timings
    python3 $S run --bruker-dda $DDA --bruker-dia $DIA --thermo $RAW \\
        --openms-bin $BUILD/bin --repeat 3 -o candidate.json
    python3 $S compare baseline.json candidate.json

    # two code paths that must give the same output (same or different report)
    python3 $S compare candidate.json \\
        --pair bruker_dia/ms2_agg+centroid@frames bruker_dia/ms2_agg+greedy2d@frames

    # tolerance comparison of two readers
    python3 $S diff --a 'thermo:run.raw:{"centroid": true}' --b mzml:run_TRFP.mzML

``compare`` exits with 0 only if every compared case is IDENTICAL, so it can
gate a commit. Hashes are taken over the exact bytes; they are meant for
before/after comparisons on one machine, not across platforms or compilers.

Requires pyopenms and numpy (``list``, ``compare`` and ``--help`` work without
them). Peak memory is not available on Windows and is reported as null there.
"""

import argparse
import datetime
import enum
import hashlib
import json
import math
import os
import platform
import re
import shutil
import signal
import socket
import statistics
import struct
import subprocess
import sys
import tempfile
import threading
import time
import traceback

try:
    import resource  # POSIX only; without it peak memory is reported as unavailable
except ImportError:
    resource = None

np = None   # numpy, imported on demand so that --help works without it
oms = None  # pyopenms, likewise

TOOL = "benchmark_vendor_readers"
REPORT_VERSION = 1
# Bump whenever the hashed byte stream changes: reports with different versions are not comparable.
FINGERPRINT_VERSION = 1
RESULT_MARKER = "@@%s-result@@ " % TOOL
ENVIRONMENT_VARIABLES = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENMS_DATA_PATH", "OPENMS_THERMO_MANAGED_DIR",
                         "OPENMS_BRUKER_SDK_PATH", "DOTNET_ROOT", "PYTHONPATH")


# ======================================================================================
# 1. The case matrix
# ======================================================================================
#
# A case is a JSON-serialisable dict:
#   {"label": str, "reader": "bruker" | "thermo" | "filehandler" | "mzml" | "topp", "input": <input name>,
#    reader-specific settings, "comment": what the case exercises}
#   bruker:      "config"            BrukerTimsFile.Config field -> value (enum fields: the member NAME)
#   thermo:      "options"           ThermoRawFileOptions field -> value
#   filehandler: "peak_file_options" see PEAK_FILE_OPTIONS
#   mzml:        "peak_file_options" (optional)
#   topp:        "command"           argument list with the placeholders {bin}, {input} and {out};
#                "out_suffix"        extension of {out}, which is then read with MzMLFile
# Inputs are given on the command line: --bruker-dda, --bruker-dia, --thermo, --input NAME=PATH.

FULL, QUICK = True, False  # run a Bruker variant on the whole file too, or only on the quick frame window

# Frame-id window of the cheap "@frames" twin that every Bruker variant gets: a few seconds from the
# middle of the gradient with several MS1 frames and the MS2 frames between them (Frames table of
# analysis.tdf: MsMsType 0 = MS1, 8 = ddaPASEF, 9 = diaPASEF).
#   DDA fixture: frames 1500-1540 = 9 MS1 + 32 MS2 frames, 205 precursors, 11 of them cut by the window
#   DIA fixture: frames 1495-1540 = 6 MS1 + 40 MS2 frames (5 cycles through all 8 window groups)
QUICK_FRAMES = {
    "bruker_dda": {"frame_id_min": 1500, "frame_id_max": 1540},
    "bruker_dia": {"frame_id_min": 1495, "frame_id_max": 1540},
}

GREEDY2D = {"ms1_centroid_mz_ppm": 10.0, "ms1_centroid_im_pct": 3.0}

# Variants that apply to both acquisition types: (variant, Config fields, whole file?, code path exercised)
BRUKER_COMMON = [
    ("default", {}, FULL,
     "AUTO export mode, all defaults. MS1: one raw spectrum per frame with a per-peak 1/K0 array. DDA MS2: one "
     "spectrum per precursor from the TOF-domain pipeline with a scalar drift time. DIA MS2: one raw spectrum "
     "per frame and isolation window with a per-peak 1/K0 array"),
    ("no_ms1", {"load_ms1": False}, FULL,
     "MS2-only load (database search): MS1 frames are not read"),
    ("no_ms2", {"load_ms2": False}, FULL,
     "MS1-only load (feature detection): MS2 frames are not read, for DDA the precursor table neither"),
    ("frame", {"export_mode": "FRAME"}, FULL,
     "loadFrames_: raw 4D frames for MS1 and MS2, no precursor grouping and no window splitting"),
    ("frame_no_ms1", {"export_mode": "FRAME", "load_ms1": False}, QUICK,
     "FRAME mode restricted to MS2 frames"),
    ("frame_no_ms2", {"export_mode": "FRAME", "load_ms2": False}, QUICK,
     "FRAME mode restricted to MS1 frames"),
    ("ms1_greedy2d", dict(GREEDY2D), FULL,
     "implicit Greedy2D MS1 centroiding: ms1_centroid_algo stays OFF, both legacy tolerances > 0 "
     "(FrameCentroider; ties of equal intensities decide which apex wins)"),
    ("ms1_greedy2d_explicit", dict(GREEDY2D, ms1_centroid_algo="GREEDY2D"), QUICK,
     "explicit ms1_centroid_algo=GREEDY2D; must equal ms1_greedy2d (compare --pair)"),
    ("ms1_greedy2d_cap", dict(GREEDY2D, ms1_centroid_max_peaks=2000), QUICK,
     "Greedy2D with the per-spectrum cap reached: only the most intense centroids are kept"),
    ("ms1_hill", {"ms1_centroid_algo": "HILL_BASED"}, FULL,
     "HillBased MS1 centroiding: hills along the IM axis within a frame (PASEFHillCentroider)"),
    ("ms1_hill_gap", {"ms1_centroid_algo": "HILL_BASED", "centroid_max_scan_gap": 1,
                      "ms1_centroid_min_hill_length": 2, "centroid_valley_factor": 1.1}, QUICK,
     "HillBased MS1 with non-default knobs: gap bridging, minimum hill length 2, eager valley splitting"),
    ("ms1_agg", {"ms1_n_neighbors": 1}, FULL,
     "MS1 frame aggregation (3-frame sum in FrameAggregator), no denoising, no centroiding; the order of "
     "equal-m/z peaks comes from a hash table followed by an unstable sort"),
    ("ms1_agg_denoise", {"ms1_n_neighbors": 2, "ms1_min_support": 2, "ms1_max_rt_distance_sec": 1.5}, QUICK,
     "5-frame MS1 aggregation with 3x3 support denoising and an RT distance cap for the neighbours"),
    ("ms1_agg+greedy2d", dict(GREEDY2D, ms1_n_neighbors=1), FULL,
     "aggregated MS1 frames, then Greedy2D centroiding"),
    ("ms1_agg+hill", {"ms1_n_neighbors": 1, "ms1_centroid_algo": "HILL_BASED"}, QUICK,
     "aggregated MS1 frames, then HillBased centroiding"),
]

BRUKER_DDA = BRUKER_COMMON + [
    ("ms2_hill", {"ms2_centroid_algo": "HILL_BASED"}, FULL,
     "DDA MS2 by IM-axis hill detection instead of the TOF-domain pipeline (default minimum hill length 2)"),
    ("ms2_hill_len1", {"ms2_centroid_algo": "HILL_BASED", "ms2_centroid_min_hill_length": 1}, QUICK,
     "HillBased DDA MS2 with the minimum hill length recommended for DDA"),
    ("ms2_greedy2d", {"ms2_centroid_algo": "GREEDY2D"}, QUICK,
     "Greedy2D does not exist for DDA MS2 and falls through to the default pipeline; must equal default"),
    ("calibrate", {"calibrate": True}, QUICK,
     "calibrate=True with the default calibration strategy: m/z comes from the MzCalibration table, the OLS "
     "recalibration is skipped; must equal default"),
    ("linear", {"tims_calibration_strategy": "LINEAR"}, QUICK,
     "linear TOF -> sqrt(m/z) and scan -> 1/K0 models from GlobalMetadata instead of the calibration tables"),
    ("calibrate_linear", {"calibrate": True, "tims_calibration_strategy": "LINEAR"}, FULL,
     "OLS m/z recalibration, which only runs on the linear model: a first pass over the precursors fits the "
     "calibration, then every frame is converted with it"),
    ("prefilter+ms1_greedy2d", dict(GREEDY2D, isotopic_prefilter=True), QUICK,
     "isotopic-partner prefilter ahead of MS1 centroiding (the prefilter never touches DDA MS2, and not raw "
     "MS1 frames either)"),
    ("ms2_hill+ms1_hill_bounds", {"ms1_centroid_algo": "HILL_BASED", "ms2_centroid_algo": "HILL_BASED",
                                  "expose_hill_bounds": True}, QUICK,
     "HillBased on both levels with expose_hill_bounds: four more float data arrays per centroided MS1 spectrum"),
]

BRUKER_DIA = BRUKER_COMMON + [
    ("ms2_hill", {"ms2_centroid_algo": "HILL_BASED"}, FULL,
     "HillBased DIA MS2 centroiding per frame and window, no RT aggregation"),
    ("ms2_agg", {"dia_ms2_n_neighbors": 1}, FULL,
     "DIA MS2 RT aggregation (3-frame sum) with the default support denoising, no centroiding: "
     "FrameAggregator::finalize, where the order of equal-m/z cells depends on the hash-table history"),
    ("ms2_agg_raw", {"dia_ms2_n_neighbors": 2, "dia_ms2_min_support": 0}, QUICK,
     "5-frame DIA MS2 sum without denoising"),
    ("ms2_agg+centroid", {"dia_ms2_n_neighbors": 1, "dia_ms2_centroid": True}, FULL,
     "legacy 2D peak picking of aggregated DIA MS2 (Gaussian smoothing + local maxima, finalizeCentroided)"),
    ("ms2_agg+greedy2d", {"dia_ms2_n_neighbors": 1, "ms2_centroid_algo": "GREEDY2D"}, QUICK,
     "the same selected by ms2_centroid_algo=GREEDY2D; must equal ms2_agg+centroid"),
    ("ms2_agg+hill", {"dia_ms2_n_neighbors": 1, "ms2_centroid_algo": "HILL_BASED"}, QUICK,
     "aggregated DIA MS2, then HillBased centroiding"),
    ("ms2_centroid_noagg", {"dia_ms2_centroid": True}, QUICK,
     "Greedy2D DIA MS2 without RT aggregation is treated as off; must equal default"),
    ("ms1_agg+ms2_agg", {"ms1_n_neighbors": 1, "dia_ms2_n_neighbors": 1, "dia_ms2_centroid": True}, QUICK,
     "MS1 and MS2 aggregation together (two neighbour windows advance through the same frames)"),
    ("isotopic_prefilter", {"isotopic_prefilter": True}, FULL,
     "isotopic-partner prefilter on raw DIA MS2 (raw MS1 frames are not prefiltered)"),
    ("prefilter+agg", {"isotopic_prefilter": True, "ms1_n_neighbors": 1, "dia_ms2_n_neighbors": 1}, QUICK,
     "prefilter after MS1 and MS2 aggregation"),
    ("hill_bounds", {"ms1_centroid_algo": "HILL_BASED", "ms2_centroid_algo": "HILL_BASED",
                     "expose_hill_bounds": True}, QUICK,
     "HillBased on both levels with expose_hill_bounds: four more float data arrays per centroided spectrum"),
]

# Range filters beyond the "@frames" twins: (label suffix, Config fields, variants to combine with, what).
# The RT window 99-105 s holds 11 MS1 + 46 MS2 frames (DDA) and 6 MS1 + 51 MS2 frames (DIA).
RT_WINDOW = {"rt_min_sec": 99.0, "rt_max_sec": 105.0}
BRUKER_RANGES = {
    "bruker_dda": [
        ("@rt", RT_WINDOW, ["default", "ms1_agg", "ms1_greedy2d", "ms2_hill"],
         "RT sub-range (resolved to a frame-id range by SQL)"),
        ("@frames+rt", {"frame_id_min": 1500, "frame_id_max": 1540, "rt_min_sec": 164.0}, ["default", "ms1_agg"],
         "frame-id range intersected with an open-ended RT range (leaves frames 1519-1540)"),
        ("@ms1_frame", {"frame_id_min": 1505, "frame_id_max": 1505}, ["default", "frame"], "a single MS1 frame"),
        ("@ms2_frame", {"frame_id_min": 1501, "frame_id_max": 1501}, ["default", "frame"],
         "a single MS2 frame: every precursor is cut by the range"),
        ("@no_frames", {"frame_id_min": 100000, "frame_id_max": 100010}, ["default"],
         "range beyond the last frame: an empty experiment with metadata"),
    ],
    "bruker_dia": [
        ("@rt", RT_WINDOW, ["default", "ms1_agg", "ms1_greedy2d", "ms2_agg+centroid"],
         "RT sub-range (resolved to a frame-id range by SQL)"),
        ("@frames+rt", {"frame_id_min": 1495, "frame_id_max": 1540, "rt_min_sec": 161.5}, ["default", "ms2_agg"],
         "frame-id range intersected with an open-ended RT range (leaves frames 1520-1540)"),
        ("@ms1_frame", {"frame_id_min": 1504, "frame_id_max": 1504}, ["default", "frame"], "a single MS1 frame"),
        ("@ms2_frame", {"frame_id_min": 1500, "frame_id_max": 1500}, ["default", "frame"],
         "a single MS2 frame: one spectrum per isolation window"),
        ("@no_frames", {"frame_id_min": 100000, "frame_id_max": 100010}, ["default"],
         "range beyond the last frame: an empty experiment with metadata"),
    ],
}

THERMO = [
    ("default", {}, "scans as acquired (profile where recorded), trailers, instrument methods and SHA-1"),
    ("centroid", {"centroid": True}, "Thermo peak picking, as FileHandler and FileConverter request it"),
    ("centroid+charge", {"centroid": True, "charge_data": True}, "adds the integer 'charge array' of the centroids"),
    ("noise", {"noise_data": True}, "adds the sampled noise/baseline arrays as meta values"),
    ("all_detectors", {"all_detectors": True}, "adds UV/PDA/analog traces and PDA spectra where the file has them"),
    ("no_trailers", {"preserve_trailers": False}, "without the per-scan trailer JSON"),
    ("no_methods", {"instrument_methods": False}, "without the embedded instrument method texts"),
    ("no_checksum", {"checksum": False}, "without the SHA-1 pass over the file"),
    ("lean", {"centroid": True, "preserve_trailers": False, "instrument_methods": False, "checksum": False},
     "the cheapest configuration: centroids only"),
    ("everything", {"centroid": True, "charge_data": True, "noise_data": True, "all_detectors": True},
     "all optional exports on"),
    ("ms1_only", {"centroid": True, "ms_levels": [1]},
     "scan selection: MS1 only; the peak data of the MSn scans is not read"),
    ("ms2_only", {"centroid": True, "ms_levels": [2]},
     "scan selection: MS2 only; the parent MS1 scans are read for the precursor intensity but not returned"),
    ("rt_window", {"centroid": True, "rt_min_sec": 100.0, "rt_max_sec": 200.0},
     "scan selection: retention time range (inclusive)"),
]

# FileHandler::loadExperiment with PeakFileOptions. Today the vendor readers load everything and FileHandler
# filters afterwards; these cases are the reference for pushing the filters down into the readers, which must
# return exactly the same. Windows per input: RT [s], m/z, intensity, precursor m/z; each cuts through the data.
FILEHANDLER_WINDOWS = {
    "bruker_dda": {"rt": [100.0, 200.0], "mz": [400.0, 800.0], "intensity": [100.0, 1e12], "precursor": [500.0, 700.0]},
    "bruker_dia": {"rt": [100.0, 200.0], "mz": [400.0, 800.0], "intensity": [100.0, 1e12], "precursor": [500.0, 700.0]},
    "thermo": {"rt": [100.0, 200.0], "mz": [400.0, 800.0], "intensity": [1e5, 1e12], "precursor": [430.0, 440.0]},
}


def _filehandler_variants(window):
    return [
        ("default", {}, "no options (Thermo: centroided, see FileHandler::loadExperiment)"),
        ("ms1", {"ms_levels": [1]}, "MS level filter: MS1 only"),
        ("ms2", {"ms_levels": [2]}, "MS level filter: MS2 only (what the search engines request)"),
        ("rt_range", {"rt_range": window["rt"]}, "RT range"),
        ("mz_range", {"mz_range": window["mz"]}, "per-peak m/z range (data arrays must be cut along)"),
        ("intensity_range", {"intensity_range": window["intensity"]}, "per-peak intensity range"),
        ("precursor_mz_range", {"precursor_mz_range": window["precursor"]}, "precursor m/z range (MS1 is kept)"),
        ("metadata_only", {"metadata_only": True}, "spectra without peaks and data arrays"),
        ("ms2+rt+precursor", {"ms_levels": [2], "rt_range": window["rt"], "precursor_mz_range": window["precursor"]},
         "several filters at once"),
    ]


# Only cheap on Thermo, so only there: options the vendor paths ignore today.
FILEHANDLER_THERMO_EXTRA = [
    ("fill_data_false", {"fill_data": False}, "fill_data=False (ignored by the vendor paths today)"),
    ("skip_chromatograms", {"skip_chromatograms": True}, "skip_chromatograms (ignored by the vendor paths today)"),
]

# End-to-end: a TOPP tool as a child process, then the fingerprint of the mzML it wrote.
TOPP = [
    ("topp/FileConverter/bruker_dda", "bruker_dda", ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}"],
     "FileConverter .d -> mzML, default parameters"),
    ("topp/FileConverter/bruker_dda_no_ms1", "bruker_dda",
     ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}", "-bruker:load_ms1", "false"],
     "FileConverter .d -> mzML, MS2 only"),
    ("topp/FileConverter/bruker_dia", "bruker_dia", ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}"],
     "FileConverter .d -> mzML, default parameters"),
    # BrukerTimsFile::transform(): the spectra are written as they are read. Each must hold the spectra of its
    # in-memory twin above (compare --pair; the two differ in what the mzML writers add around the spectra)
    ("topp/FileConverter/bruker_dda_lowmemory", "bruker_dda",
     ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}", "-process_lowmemory"],
     "FileConverter .d -> mzML without loading the run into memory"),
    ("topp/FileConverter/bruker_dda_no_ms1_lowmemory", "bruker_dda",
     ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}", "-process_lowmemory", "-bruker:load_ms1", "false"],
     "FileConverter .d -> mzML without loading the run into memory, MS2 only"),
    ("topp/FileConverter/bruker_dia_lowmemory", "bruker_dia",
     ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}", "-process_lowmemory"],
     "FileConverter .d -> mzML without loading the run into memory"),
    ("topp/FileConverter/thermo", "thermo",
     ["{bin}/FileConverter", "-in", "{input}", "-out", "{out}", "-RawToMzML:reader", "inprocess"],
     "FileConverter .raw -> mzML with the in-process reader (vendor peak picking)"),
]


def _default_cases():
    cases = []
    for input_name, variants in (("bruker_dda", BRUKER_DDA), ("bruker_dia", BRUKER_DIA)):
        frames = QUICK_FRAMES[input_name]
        window = "frames %d-%d" % (frames["frame_id_min"], frames["frame_id_max"])
        configs = {}
        for variant, config, whole_file, what in variants:
            configs[variant] = config
            label = "%s/%s" % (input_name, variant)
            if whole_file:
                cases.append({"label": label, "reader": "bruker", "input": input_name, "config": dict(config),
                              "comment": what})
            cases.append({"label": label + "@frames", "reader": "bruker", "input": input_name,
                          "config": dict(config, **frames), "comment": "%s [%s]" % (what, window)})
        for suffix, range_config, combine_with, what in BRUKER_RANGES[input_name]:
            for variant in combine_with:
                cases.append({"label": "%s/%s%s" % (input_name, variant, suffix), "reader": "bruker",
                              "input": input_name, "config": dict(configs[variant], **range_config),
                              "comment": "%s; variant %s" % (what, variant)})
    for variant, options, what in THERMO:
        cases.append({"label": "thermo/" + variant, "reader": "thermo", "input": "thermo", "options": dict(options),
                      "comment": what})
    for input_name in ("bruker_dda", "bruker_dia", "thermo"):
        variants = _filehandler_variants(FILEHANDLER_WINDOWS[input_name])
        if input_name == "thermo":
            variants += FILEHANDLER_THERMO_EXTRA
        for variant, options, what in variants:
            cases.append({"label": "filehandler/%s/%s" % (input_name, variant), "reader": "filehandler",
                          "input": input_name, "peak_file_options": dict(options), "comment": what})
    for label, input_name, command, what in TOPP:
        cases.append({"label": label, "reader": "topp", "input": input_name, "command": list(command),
                      "out_suffix": ".mzML", "comment": what})
    cases.append({"label": "mzml/reference", "reader": "mzml", "input": "reference_mzml",
                  "comment": "fingerprint of an existing mzML (--input reference_mzml=FILE)"})
    return cases


DEFAULT_CASES = _default_cases()
READERS = ("bruker", "thermo", "filehandler", "mzml", "topp")
SETTING_KEYS = ("config", "options", "peak_file_options", "command", "out_suffix")

# PeakFileOptions keys understood in "peak_file_options": key -> (setter, kind of value)
PEAK_FILE_OPTIONS = {
    "ms_levels": ("setMSLevels", "ints"),
    "rt_range": ("setRTRange", "range"),
    "mz_range": ("setMZRange", "range"),
    "intensity_range": ("setIntensityRange", "range"),
    "precursor_mz_range": ("setPrecursorMZRange", "range"),
    "metadata_only": ("setMetadataOnly", "bool"),
    "fill_data": ("setFillData", "bool"),
    "skip_chromatograms": ("setSkipChromatograms", "bool"),
    "sort_spectra_by_mz": ("setSortSpectraByMZ", "bool"),
    "precursor_mz_selected_ion": ("setPrecursorMZSelectedIon", "bool"),
}


class CaseSkipped(Exception):
    """The case cannot run here (input not given, feature not in this build); not a failure."""


class CaseError(Exception):
    """The case is wrong (unknown enum member, bad value) or the reader failed."""


def load_cases(path):
    """The case matrix: the built-in one, or the JSON file `path` (a list of cases, or {"cases": [...]})."""
    if path is None:
        cases = [dict(case) for case in DEFAULT_CASES]
    else:
        with open(path) as handle:
            cases = json.load(handle)
        if isinstance(cases, dict):
            cases = cases.get("cases")
        if not isinstance(cases, list):
            raise SystemExit("%s: expected a JSON list of cases" % path)
    seen = set()
    for case in cases:
        if not isinstance(case, dict) or not isinstance(case.get("label"), str):
            raise SystemExit("case without a 'label': %r" % (case,))
        label = case["label"]
        if label in seen:
            raise SystemExit("duplicate case label '%s'" % label)
        seen.add(label)
        if case.get("reader") not in READERS:
            raise SystemExit("case '%s': 'reader' must be one of %s" % (label, ", ".join(READERS)))
        if not isinstance(case.get("input"), str):
            raise SystemExit("case '%s': 'input' must name an input" % label)
        if case["reader"] == "topp" and not isinstance(case.get("command"), list):
            raise SystemExit("case '%s': a topp case needs a 'command' argument list" % label)
    return cases


def select_cases(cases, include, exclude):
    try:
        if include:
            cases = [case for case in cases if re.search(include, case["label"])]
        if exclude:
            cases = [case for case in cases if not re.search(exclude, case["label"])]
    except re.error as error:
        raise SystemExit("invalid regular expression: %s" % error)
    return cases


def case_settings(case):
    """The reader-specific part of a case (what distinguishes it besides reader and input)."""
    return {key: case[key] for key in SETTING_KEYS if key in case}


# ======================================================================================
# 2. pyOpenMS access: what the fingerprint reads, and what this build provides
# ======================================================================================

def _import_runtime_deps():
    global np, oms
    if oms is not None:
        return
    try:
        import numpy as _np
    except ImportError as error:
        raise SystemExit("missing dependency: %s (pip install numpy)" % error)
    try:
        import pyopenms as _oms
    except ImportError as error:
        raise SystemExit("cannot import pyopenms (%s): build OpenMS-build/pyOpenMS and set "
                         "PYTHONPATH=OpenMS-build/pyOpenMS, or `pip install pyopenms`" % error)
    np, oms = _np, _oms


# Everything the fingerprint reads from the metadata objects, by pyOpenMS accessor name. The table is the
# specification of the metadata digests:
#   values    hashed exactly (floats by their bytes)
#   paths     depend on where the input lies; collected in "path_dependent" and hashed into the separate
#             component "source_paths", so that the same data reached through another path or a symlink
#             differs in that one component only (compare --ignore-paths)
#   volatile  differ from build to build or run to run (software version with build date, completion time);
#             recorded under "volatile" for information and never part of a verdict
#   children  accessor -> entry of this table describing the returned object, list of objects or dict of objects
#   items     entry describing the elements of an indexable container (size() / [i])
#   meta      True: all meta values (sorted key, type-tagged exact value); "volatile": recorded, not hashed
#   class     pyOpenMS class if it is not the entry name
# An accessor that the installed pyOpenMS does not have is detected once at start-up (detect_coverage), listed
# under "fingerprint_coverage"/"missing" in the report header, and left out; nothing is skipped silently.
SCHEMA = {
    "MSSpectrum": {
        "values": ["getMSLevel", "getRT", "getNativeID", "getName", "getComment", "getType", "getDriftTime",
                   "getDriftTimeUnit", "getIMFormat", "getIMPeakType"],
        "children": [("getInstrumentSettings", "InstrumentSettings"), ("getAcquisitionInfo", "AcquisitionInfo"),
                     ("getPrecursors", "Precursor"), ("getProducts", "Product"), ("getSourceFile", "SourceFile"),
                     ("getDataProcessing", "DataProcessing")],
        "meta": True,
    },
    "MSChromatogram": {
        "values": ["getNativeID", "getName", "getComment", "getChromatogramType", "getMZ"],
        "children": [("getPrecursor", "Precursor"), ("getProduct", "Product"),
                     ("getInstrumentSettings", "InstrumentSettings"), ("getAcquisitionInfo", "AcquisitionInfo"),
                     ("getSourceFile", "SourceFile"), ("getDataProcessing", "DataProcessing")],
        "meta": True,
    },
    "InstrumentSettings": {
        "values": ["getScanMode", "getZoomScan", "getPolarity"],
        "children": [("getScanWindows", "ScanWindow")],
        "meta": True,
    },
    "ScanWindow": {"values": ["begin", "end"], "meta": True},
    "AcquisitionInfo": {"values": ["getMethodOfCombination"], "items": "Acquisition", "meta": True},
    "Acquisition": {"values": ["getIdentifier"], "meta": True},
    "Precursor": {
        "values": ["getMZ", "getCharge", "getIntensity", "getIsolationWindowLowerOffset",
                   "getIsolationWindowUpperOffset", "getDriftTime", "getDriftTimeUnit",
                   "getDriftTimeWindowLowerOffset", "getDriftTimeWindowUpperOffset", "getActivationMethods",
                   "getActivationEnergy", "getPossibleChargeStates", "getCVTerms"],
        "meta": True,
    },
    "Product": {
        "values": ["getMZ", "getIsolationWindowLowerOffset", "getIsolationWindowUpperOffset", "getCVTerms"],
        "meta": True,
    },
    "SourceFile": {
        "values": ["getFileType", "getFileSize", "getChecksum", "getChecksumType", "getNativeIDType",
                   "getNativeIDTypeAccession", "getCVTerms"],
        "paths": ["getNameOfFile", "getPathToFile"],
        "meta": True,
    },
    "DataProcessing": {
        "values": ["getProcessingActions"],
        "volatile": ["getCompletionTime"],
        "children": [("getSoftware", "ProcessingSoftware")],
        "meta": "volatile",  # TOPP tools store their parameters here, input and output paths included
    },
    "ProcessingSoftware": {"class": "Software", "values": ["getName"], "volatile": ["getVersion"]},
    "MSExperiment": {
        "values": ["getComment", "getFractionIdentifier", "getIdentifier", "getLoadedFileType", "getDateTime"],
        "paths": ["getLoadedFilePath", "getPrimaryMSRunPath"],
        "children": [("getSourceFiles", "SourceFile"), ("getInstrument", "Instrument"),
                     ("getInstrumentConfigurations", "Instrument"), ("getSample", "Sample"),
                     ("getContacts", "ContactPerson"), ("getHPLC", "HPLC")],
        "meta": True,
    },
    "Instrument": {
        "values": ["getName", "getVendor", "getModel", "getCustomizations", "getIonOptics"],
        "children": [("getSoftware", "Software"), ("getIonSources", "IonSource"),
                     ("getMassAnalyzers", "MassAnalyzer"), ("getIonDetectors", "IonDetector")],
        "meta": True,
    },
    "Software": {"values": ["getName", "getVersion"]},  # acquisition software, read from the file
    "IonSource": {"values": ["getInletType", "getIonizationMethod", "getPolarity", "getOrder"], "meta": True},
    "MassAnalyzer": {
        "values": ["getType", "getResolutionMethod", "getResolutionType", "getScanDirection", "getScanLaw",
                   "getReflectronState", "getResolution", "getAccuracy", "getScanRate", "getScanTime",
                   "getTOFTotalPathLength", "getIsolationWidth", "getFinalMSExponent", "getMagneticFieldStrength",
                   "getOrder"],
        "meta": True,
    },
    "IonDetector": {
        "values": ["getType", "getAcquisitionMode", "getResolution", "getADCSamplingFrequency", "getOrder"],
        "meta": True,
    },
    "Sample": {
        "values": ["getName", "getOrganism", "getNumber", "getComment", "getState", "getMass", "getVolume",
                   "getConcentration"],
        "children": [("getSubsamples", "Sample")],
        "meta": True,
    },
    "ContactPerson": {
        "values": ["getFirstName", "getLastName", "getInstitution", "getEmail", "getURL", "getAddress",
                   "getContactInfo"],
        "meta": True,
    },
    "HPLC": {
        "values": ["getInstrument", "getColumn", "getTemperature", "getPressure", "getFlux", "getComment"],
        "children": [("getGradient", "Gradient")],
    },
    "Gradient": {"values": ["getEluents", "getTimepoints", "getPercentages"]},
}

# Value ranges kept by the experiment (updateRanges()); an empty range raises and is hashed as None.
EXPERIMENT_RANGES = ("getMinRT", "getMaxRT", "getMinMZ", "getMaxMZ", "getMinIntensity", "getMaxIntensity",
                     "getMinMobility", "getMaxMobility")

# (component, list accessor on MSSpectrum / MSChromatogram, element class)
DATA_ARRAY_KINDS = (
    ("float_arrays", "getFloatDataArrays", "FloatDataArray"),
    ("integer_arrays", "getIntegerDataArrays", "IntegerDataArray"),
    ("string_arrays", "getStringDataArrays", "StringDataArray"),
)

# Known gaps that no pyOpenMS accessor closes; printed with the coverage so nobody assumes they are hashed.
NOT_REACHABLE = [
    "unit of a meta value (DataValue unit) and unit of a CV term",
    "data processing attached to a single float/integer/string data array",
    "spectrum-level identifications (PeptideIdentification) and experiment-level protein identifications",
]

_coverage_cache = None


def detect_coverage(pinned=None):
    """Work out ONCE which accessors of SCHEMA this pyOpenMS build provides.

    Returns {"use": what the fingerprint functions read, "report": the "fingerprint_coverage" header entry}.
    `pinned` is the "hashed" table of an earlier report (run --coverage-from): accessors that are not in it
    are left out too, so that a build with more bindings still yields digests comparable with that report.
    """
    global _coverage_cache
    if _coverage_cache is not None and pinned is None:
        return _coverage_cache
    _import_runtime_deps()
    use, hashed, missing, unpinned = {}, {}, {}, {}

    def consider(entry_name, description, available):
        if not available:
            missing.setdefault(entry_name, []).append(description)
        elif pinned is not None and description not in pinned.get(entry_name, []):
            unpinned.setdefault(entry_name, []).append(description)
        else:
            hashed.setdefault(entry_name, []).append(description)
            return True
        return False

    for entry_name, spec in SCHEMA.items():
        cls = getattr(oms, spec.get("class", entry_name), None)
        entry = {"values": [], "paths": [], "volatile": [], "children": [], "items": None, "meta": False}
        for kind in ("values", "paths", "volatile"):
            for accessor in spec.get(kind, []):
                description = accessor if kind == "values" else "%s (%s)" % (accessor, kind)
                if consider(entry_name, description, cls is not None and hasattr(cls, accessor)):
                    entry[kind].append(accessor)
        for accessor, child in spec.get("children", []):
            if consider(entry_name, "%s -> %s" % (accessor, child), cls is not None and hasattr(cls, accessor)):
                entry["children"].append((accessor, child))
        if "items" in spec:
            available = cls is not None and hasattr(cls, "size") and hasattr(cls, "__getitem__")
            if consider(entry_name, "[i] -> %s" % spec["items"], available):
                entry["items"] = spec["items"]
        if spec.get("meta"):
            description = "meta values" if spec["meta"] is True else "meta values (volatile)"
            available = cls is not None and hasattr(cls, "getKeys") and hasattr(cls, "getMetaValue")
            if consider(entry_name, description, available):
                entry["meta"] = spec["meta"]
        use[entry_name] = entry

    experiment = getattr(oms, "MSExperiment", None)
    use["ranges"] = [name for name in EXPERIMENT_RANGES
                     if consider("MSExperiment", "%s (range)" % name, hasattr(experiment, name))]
    use["arrays"] = {}
    for component, accessor, element in DATA_ARRAY_KINDS:
        element_class = getattr(oms, element, None)
        available = (all(hasattr(getattr(oms, owner, None), accessor) for owner in ("MSSpectrum", "MSChromatogram"))
                     and all(hasattr(element_class, name) for name in ("getName", "get_data")))
        use["arrays"][component] = consider(
            element, "MSSpectrum/MSChromatogram.%s: getName, get_data" % accessor, available)
        array_meta = hasattr(element_class, "getKeys") and hasattr(element_class, "getMetaValue")
        use["arrays"][component + "_meta"] = use["arrays"][component] and consider(element, "meta values", array_meta)
    for owner in ("MSSpectrum", "MSChromatogram"):
        if not hasattr(getattr(oms, owner, None), "get_peaks"):
            raise SystemExit("this pyOpenMS has no %s.get_peaks(): no fingerprint possible" % owner)
        hashed.setdefault(owner, []).append("get_peaks (float64 positions, float32 intensities)")
    use["spectrum_view"] = hasattr(experiment, "spectrum_view")

    report = {"fingerprint_version": FINGERPRINT_VERSION, "hashed": hashed, "missing": missing,
              "not_reachable_from_pyopenms": NOT_REACHABLE,
              "spectrum_access": "MSExperiment.spectrum_view (no copy)" if use["spectrum_view"]
                                 else "MSExperiment[i] (copy)"}
    if pinned is not None:
        report["left_out_to_match_pinned_coverage"] = unpinned
    coverage = {"use": use, "report": report}
    if pinned is None:
        _coverage_cache = coverage
    return coverage


def _get(obj, name):
    value = getattr(obj, name)
    return value() if callable(value) else value


def _plain(value):
    """Turn the few pyOpenMS objects that accessors return into plain data; everything else is left alone."""
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    kind = type(value).__name__
    if kind == "DateTime":
        return value.get()
    if kind == "CVTerm":
        return [value.getAccession(), value.getName(), value.getCVIdentifierRef(),
                value.getValue() if value.hasValue() else None]
    return value


def _meta_items(obj):
    """All meta values of an object as [[key, value], ...], sorted by key."""
    try:
        keys = obj.getKeys()
    except TypeError:  # older signature: fills the list it is given
        keys = []
        obj.getKeys(keys)
    return [[key, obj.getMetaValue(key)] for key in sorted(keys)]


def _note(counter, label, value):
    """Count a value that is recorded for information (path-dependent or volatile), not hashed with the rest."""
    key = (label, json.dumps(_jsonable(value), sort_keys=True))
    counter[key] = counter.get(key, 0) + 1


def _jsonable(value):
    """A JSON-serialisable rendering for the report (the digests never use this)."""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = [_jsonable(item) for item in value]
        return sorted(items, key=repr) if isinstance(value, (set, frozenset)) else items
    if isinstance(value, enum.Enum):
        return value.name
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)  # 'inf', 'nan': keeps the report strict JSON
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _describe(obj, entry_name, use, side):
    """Walk one metadata object as laid out in SCHEMA.

    Returns (values, meta_values): two nested lists for _encode(). Path-dependent and volatile values go to
    the counters side["paths"] and side["volatile"] instead.
    """
    entry = use[entry_name]
    values = [[name, _plain(_get(obj, name))] for name in entry["values"]]
    for name in entry["paths"]:
        value = _get(obj, name)
        if value:  # unset names and paths are left out: the per-spectrum source file is normally empty
            _note(side["paths"], "%s.%s" % (entry_name, name), value)
    for name in entry["volatile"]:
        _note(side["volatile"], "%s.%s" % (entry_name, name), _plain(_get(obj, name)))
    meta_values = []
    if entry["meta"] is True:
        meta_values.append(["", _meta_items(obj)])
    elif entry["meta"] == "volatile":
        _note(side["volatile"], "%s meta values" % entry_name, _meta_items(obj))
    groups = []
    for accessor, child_name in entry["children"]:
        child = _get(obj, accessor)
        if isinstance(child, dict):
            members = sorted(child.items())
        elif isinstance(child, (list, tuple)):
            members = [(None, item) for item in child]
        else:
            members = [(None, child)]
        groups.append((accessor, child_name, members))
    if entry["items"]:
        groups.append(("[i]", entry["items"], [(None, obj[i]) for i in range(obj.size())]))
    for accessor, child_name, members in groups:
        child_values, child_meta = [], []
        for key, member in members:
            member_values, member_meta = _describe(member, child_name, use, side)
            child_values.append([key, member_values])
            child_meta.append([key, member_meta])
        values.append([accessor, child_values])
        meta_values.append([accessor, child_meta])
    return values, meta_values


# ======================================================================================
# 3. The fingerprint
# ======================================================================================

def _encode(value):
    """Exact, unambiguous, type-tagged byte encoding of a (nested) metadata value.

    Floats are written as their 8 IEEE-754 bytes, never through repr(); strings and containers carry their
    length; the tag keeps 1, 1.0, "1" and True apart. Anything without an exact encoding raises: a value the
    fingerprint cannot represent must not go unnoticed.
    """
    if value is None:
        return b"N"
    if value is True:
        return b"T"
    if value is False:
        return b"F"
    kind = type(value)
    if kind is int:
        text = b"%d" % value
        return b"i%d:%s" % (len(text), text)
    if kind is float:
        return b"d" + struct.pack("<d", value)
    if kind is str:
        data = value.encode("utf-8", "surrogatepass")
        return b"s%d:%s" % (len(data), data)
    if kind is bytes:
        return b"b%d:%s" % (len(value), value)
    if kind in (list, tuple):
        if value and all(type(item) is float for item in value):  # long float lists (sampled noise arrays)
            return b"D%d:%s" % (len(value), struct.pack("<%dd" % len(value), *value))
        return b"l%d:%s" % (len(value), b"".join(_encode(item) for item in value))
    if kind is dict:
        items = sorted((_encode(key), _encode(item)) for key, item in value.items())
        return b"m%d:%s" % (len(items), b"".join(key + item for key, item in items))
    if kind in (set, frozenset):
        items = sorted(_encode(item) for item in value)
        return b"S%d:%s" % (len(items), b"".join(items))
    if isinstance(value, enum.Enum):  # pyOpenMS enums: by class and member name, not by number
        return b"e" + _encode("%s.%s" % (kind.__qualname__, value.name))
    if np is not None and isinstance(value, (np.ndarray, np.generic)):
        array = np.ascontiguousarray(value)
        return b"a" + _encode(array.dtype.str) + _encode(list(array.shape)) + array.tobytes()
    raise TypeError("fingerprint: no exact encoding for a value of type %s.%s" % (kind.__module__, kind.__qualname__))


def _digest(data):
    return hashlib.sha256(data).digest()


def _array_bytes(array):
    """dtype, length and raw bytes of a numeric array: what is hashed for every peak or data array."""
    return b"%s:%d:%s" % (array.dtype.str.encode(), len(array), array.tobytes())


def _read_data_arrays(obj, use):
    """The data arrays of a spectrum or chromatogram: {component: [(name, meta values, data), ...] or None}.

    Numeric data comes as a contiguous numpy array (float32 / int32), string data as a list of str. None means
    this pyOpenMS build cannot read that kind of array (see the coverage in the report header).
    """
    arrays = {}
    for component, accessor, _ in DATA_ARRAY_KINDS:
        if not use["arrays"][component]:
            arrays[component] = None
            continue
        entries = []
        for array in getattr(obj, accessor)():
            data = array.get_data()
            if component != "string_arrays":
                data = np.ascontiguousarray(data)
            meta = _meta_items(array) if use["arrays"][component + "_meta"] else None
            entries.append((array.getName(), meta, data))
        arrays[component] = entries
    return arrays


def _data_array_digests(arrays, order, n_points):
    """One digest per kind of data array. With `order`, arrays as long as the peak list are permuted along."""
    digests = {}
    for component, entries in arrays.items():
        if entries is None:
            digests[component] = None
            continue
        stream = hashlib.sha256(b"%d:" % len(entries))  # the number of arrays is part of the stream
        for name, meta, data in entries:
            stream.update(_encode([name, meta]))
            permute = order is not None and len(data) == n_points
            if isinstance(data, list):
                stream.update(_encode([data[i] for i in order] if permute else data))
            else:
                stream.update(_array_bytes(data[order] if permute else data))
        digests[component] = stream.digest()
    return digests


def _peak_digests(positions, intensities, arrays, order=None):
    """Digests of the peak data of one spectrum, in stored order or permuted by `order`."""
    digests = _data_array_digests(arrays, order, len(positions))
    if order is not None:
        positions, intensities = positions[order], intensities[order]
    digests["mz"] = _digest(_array_bytes(positions))
    digests["intensity"] = _digest(_array_bytes(intensities))
    return digests


def _bits(array):
    """The values as unsigned integers with the same bytes.

    Sorting by bit pattern gives a total order in which every distinct value has its own place (0.0 and -0.0,
    NaN payloads), and it is the numeric order for positive floats.
    """
    if array.dtype.kind in "fiu" and array.dtype.itemsize in (1, 2, 4, 8):
        return array.view("u%d" % array.dtype.itemsize)
    return array


def _canonical_order(positions, intensities, arrays):
    """Permutation into the canonical peak order, or None if the peaks already are in it.

    Canonical order: by m/z, then by every float data array (ion mobility), then by intensity, then by the
    integer and string data arrays; only data arrays with one value per peak take part. Peaks that still tie
    are equal in everything that is hashed, so their order does not matter.
    """
    n = len(positions)
    if n < 2:
        return None
    keys = [_bits(positions)]
    keys += [_bits(data) for _, _, data in arrays["float_arrays"] or [] if len(data) == n]
    keys.append(_bits(intensities))
    keys += [_bits(data) for _, _, data in arrays["integer_arrays"] or [] if len(data) == n]
    keys += [np.unique(np.array(data), return_inverse=True)[1] for _, _, data in arrays["string_arrays"] or []
             if len(data) == n]
    tied = np.ones(n - 1, dtype=bool)  # neighbours that are equal in all keys looked at so far
    for key in keys:
        if np.any(tied & (key[:-1] > key[1:])):
            order = _order_within_equal_mz(keys)
            return np.lexsort(keys[::-1]) if order is None else order  # lexsort takes the primary key last
        tied &= key[:-1] == key[1:]
        if not tied.any():
            break
    return None


def _order_within_equal_mz(keys):
    """Shortcut for _canonical_order() that covers what the Bruker reader returns, or None if it does not apply.

    There the peaks are sorted by m/z already and only peaks of equal m/z (the same TOF bin in several mobility
    scans, more than half of an MS1 frame) have to be put in order by the next key. Sorting one 64-bit number,
    (rank of the m/z value, second key), is then enough and, being nearly sorted, an order of magnitude faster
    than the general sort over all keys. The result is the same permutation.
    """
    first, second = keys[0], keys[1]
    if second.dtype.kind != "u" or second.dtype.itemsize > 4 or np.any(first[:-1] > first[1:]):
        return None
    rank = np.zeros(len(first), dtype=np.uint64)
    np.cumsum(first[1:] != first[:-1], out=rank[1:])
    combined = (rank << np.uint64(32)) | second
    order = np.argsort(combined, kind="stable")
    combined = combined[order]
    if len(keys) > 2 and np.any(combined[1:] == combined[:-1]):
        return None  # peaks equal in both keys: the remaining keys decide
    return order


SPECTRUM_COMPONENTS = ("mz", "intensity", "float_arrays", "integer_arrays", "string_arrays", "spectrum_meta",
                       "spectrum_meta_values")
# Every digest that must be equal for the verdict IDENTICAL, in the order they are reported.
ORDER_SENSITIVE = SPECTRUM_COMPONENTS + ("peak_counts", "chromatograms", "settings", "ranges", "source_paths")


def _combine(parts):
    """One digest for a spectrum from its component digests (a component that cannot be read counts as empty)."""
    return _digest(b"".join(parts[name] or b"-" for name in SPECTRUM_COMPONENTS))


def _count_arrays(statistics_by_name, arrays):
    for component, entries in arrays.items():
        for name, _, data in entries or []:
            entry = statistics_by_name.setdefault(component, {}).setdefault(name, {"arrays": 0, "values": 0})
            entry["arrays"] += 1
            entry["values"] += len(data)


def _replace_paths(value, aliases):
    if isinstance(value, str):
        for real, alias in aliases:
            value = value.replace(real, alias)
        return value
    if isinstance(value, list):
        return [_replace_paths(item, aliases) for item in value]
    return value


def _noted(counter, aliases=()):
    """The values counted by _note() as a sorted list of [label, value, occurrences]."""
    merged = {}
    for (label, text), count in counter.items():
        key = (label, json.dumps(_replace_paths(json.loads(text), aliases), sort_keys=True))
        merged[key] = merged.get(key, 0) + count
    return [[label, json.loads(text), count] for (label, text), count in sorted(merged.items())]


def fingerprint_experiment(exp, per_spectrum=True, path_aliases=None, coverage=None):
    """Content fingerprint of an in-memory MSExperiment.

    Returns a JSON-serialisable dict:
      counts          spectra, peaks (total and per MS level), chromatograms, data arrays by name, precursors ...
      digests         one SHA-256 per component over all spectra in container order; all of them equal means
                      the output is IDENTICAL: mz, intensity, float_arrays (ion mobility lives here),
                      integer_arrays, string_arrays, spectrum_meta, spectrum_meta_values, peak_counts,
                      chromatograms, settings, ranges, source_paths
      canonical       digests that ignore the order of peaks within a spectrum and/or the order of spectra,
                      to tell a re-ordering from a change of content
      per_spectrum    [native ID, MS level, RT, peak count, 16 hex digits of the digest of everything hashed
                      for the spectrum, 16 hex digits of the same with the peaks in canonical order]
      path_dependent  file names and paths found in the metadata (hashed into source_paths only)
      volatile        software versions, completion times: recorded, never part of a verdict

    Each component digest is a hash over the per-spectrum digests of that component, and each of those covers
    the array length, so a peak that moves to the neighbouring spectrum changes the digest although the
    concatenated arrays stay the same. Spectra are processed one at a time; memory stays bounded by the
    largest spectrum.

    `path_aliases` maps path prefixes (temporary directories) to placeholders in the path-dependent values.
    """
    _import_runtime_deps()
    use = (coverage or detect_coverage())["use"]
    aliases = sorted((path_aliases or {}).items(), key=lambda item: -len(item[0]))
    side = {"paths": {}, "volatile": {}}
    streams = {name: hashlib.sha256() for name in SPECTRUM_COMPONENTS}
    peak_count_stream = hashlib.sha256()
    ordered, canonical, records = [], [], []
    counts = {"spectra": 0, "peaks": 0, "spectra_by_ms_level": {}, "peaks_by_ms_level": {}, "empty_spectra": 0,
              "precursors": 0, "products": 0, "spectra_with_scalar_drift_time": 0, "spectra_not_in_canonical_peak_order": 0}
    array_statistics = {}

    n_spectra = exp.getNrSpectra()
    get_spectrum = exp.spectrum_view if use["spectrum_view"] else exp.__getitem__
    for index in range(n_spectra):
        spectrum = get_spectrum(index)
        positions, intensities = spectrum.get_peaks()  # float64 m/z, float32 intensity (copies)
        positions, intensities = np.ascontiguousarray(positions), np.ascontiguousarray(intensities)
        arrays = _read_data_arrays(spectrum, use)
        values, meta_values = _describe(spectrum, "MSSpectrum", use, side)
        fields = dict(values)

        parts = _peak_digests(positions, intensities, arrays)
        parts["spectrum_meta"] = _digest(_encode(values))
        parts["spectrum_meta_values"] = _digest(_encode(meta_values))
        for name in SPECTRUM_COMPONENTS:
            streams[name].update(parts[name] or b"-")
        peak_count_stream.update(struct.pack("<Q", len(positions)))
        ordered_digest = _combine(parts)

        order = _canonical_order(positions, intensities, arrays)
        if order is None:
            canonical_digest = ordered_digest
        else:
            counts["spectra_not_in_canonical_peak_order"] += 1
            canonical_digest = _combine(dict(parts, **_peak_digests(positions, intensities, arrays, order)))
        ordered.append(ordered_digest)
        canonical.append(canonical_digest)

        level = str(fields["getMSLevel"])
        counts["spectra"] += 1
        counts["peaks"] += len(positions)
        counts["spectra_by_ms_level"][level] = counts["spectra_by_ms_level"].get(level, 0) + 1
        counts["peaks_by_ms_level"][level] = counts["peaks_by_ms_level"].get(level, 0) + len(positions)
        counts["empty_spectra"] += len(positions) == 0
        counts["precursors"] += len(fields.get("getPrecursors", []))
        counts["products"] += len(fields.get("getProducts", []))
        counts["spectra_with_scalar_drift_time"] += fields.get("getDriftTime", -1.0) != -1.0
        _count_arrays(array_statistics, arrays)
        if per_spectrum:
            records.append([fields["getNativeID"], fields["getMSLevel"], _jsonable(fields["getRT"]), len(positions),
                            ordered_digest.hex()[:16], canonical_digest.hex()[:16]])

    chromatogram_digests = []
    counts["chromatograms"] = exp.getNrChromatograms()
    counts["chromatogram_points"] = 0
    for index in range(counts["chromatograms"]):
        chromatogram = exp.getChromatogram(index)
        times, intensities = chromatogram.get_peaks()  # float64 RT, float32 intensity
        times, intensities = np.ascontiguousarray(times), np.ascontiguousarray(intensities)
        arrays = _read_data_arrays(chromatogram, use)
        values, meta_values = _describe(chromatogram, "MSChromatogram", use, side)
        array_digests = _data_array_digests(arrays, None, len(times))
        chromatogram_digests.append(_digest(
            _encode([values, meta_values]) + _array_bytes(times) + _array_bytes(intensities)
            + b"".join(array_digests[component] or b"-" for component, _, _ in DATA_ARRAY_KINDS)))
        counts["chromatogram_points"] += len(times)
        _count_arrays(array_statistics, arrays)
    counts["data_arrays"] = array_statistics

    settings_values, settings_meta = _describe(exp, "MSExperiment", use, side)
    ranges = []
    for name in use["ranges"]:
        try:
            ranges.append([name, getattr(exp, name)()])
        except RuntimeError:  # empty range
            ranges.append([name, None])

    path_dependent = _noted(side["paths"], aliases)
    digests = {name: stream.hexdigest() for name, stream in streams.items()}
    for component, _, _ in DATA_ARRAY_KINDS:
        if not use["arrays"][component]:
            digests[component] = None  # not readable with this pyOpenMS build, see the coverage
    digests["peak_counts"] = peak_count_stream.hexdigest()
    digests["chromatograms"] = _digest(b"".join(chromatogram_digests)).hex()
    digests["settings"] = _digest(_encode([settings_values, settings_meta])).hex()
    digests["ranges"] = _digest(_encode(ranges)).hex()
    digests["source_paths"] = _digest(_encode(path_dependent)).hex()
    fingerprint = {
        "fingerprint_version": FINGERPRINT_VERSION,
        "counts": counts,
        "digests": digests,
        "canonical": {
            # peaks of every spectrum in canonical order, spectra in container order
            "ignoring_peak_order": _digest(b"".join(canonical)).hex(),
            # peaks in stored order, spectra as a sorted multiset
            "ignoring_spectrum_order": _digest(b"".join(sorted(ordered))).hex(),
            # both: changes only if the content of some spectrum changes
            "ignoring_both": _digest(b"".join(sorted(canonical))).hex(),
            "chromatograms_ignoring_order": _digest(b"".join(sorted(chromatogram_digests))).hex(),
        },
        "ranges": _jsonable(ranges),
        "path_dependent": path_dependent,
        "volatile": _noted(side["volatile"]),
    }
    if per_spectrum:
        fingerprint["per_spectrum"] = records
    return fingerprint


# ======================================================================================
# 4. Comparing two fingerprints
# ======================================================================================

def _flatten(counts, prefix=""):
    flat = {}
    for key, value in counts.items():
        if isinstance(value, dict):
            flat.update(_flatten(value, prefix + key + "."))
        else:
            flat[prefix + key] = value
    return flat


def _record(index, record):
    return {"index": index, "native_id": record[0], "ms_level": record[1], "rt": record[2], "peaks": record[3]}


def _spectrum_differences(base, cand, limit):
    """Name the spectra that differ, from the per-spectrum records of two fingerprints."""
    if base is None or cand is None:
        return {"available": False}
    result = {"available": True, "aligned_by": "index", "n_base": len(base), "n_cand": len(cand)}
    base_ids, cand_ids = [record[0] for record in base], [record[0] for record in cand]
    if base_ids == cand_ids:
        pairs = [(i, i) for i in range(len(base))]
    elif len(set(base_ids)) == len(base_ids) and len(set(cand_ids)) == len(cand_ids):
        # different sequence of native IDs, but they are unique: pair the spectra by ID
        result["aligned_by"] = "native ID"
        position = {native_id: i for i, native_id in enumerate(cand_ids)}
        pairs = [(i, position[native_id]) for i, native_id in enumerate(base_ids) if native_id in position]
        in_base = set(base_ids)
        only_base = [i for i, native_id in enumerate(base_ids) if native_id not in position]
        only_cand = [i for i, native_id in enumerate(cand_ids) if native_id not in in_base]
        result["n_only_in_base"], result["n_only_in_cand"] = len(only_base), len(only_cand)
        result["only_in_base"] = [_record(i, base[i]) for i in only_base[:limit]]
        result["only_in_cand"] = [_record(i, cand[i]) for i in only_cand[:limit]]
        moved = [(i, j) for (i, j), (_, previous) in zip(pairs[1:], pairs[:-1]) if j < previous]
        result["n_out_of_order"] = len(moved)
        result["out_of_order"] = [dict(_record(i, base[i]), index_cand=j) for i, j in moved[:limit]]
    elif len(base) == len(cand):
        pairs = [(i, i) for i in range(len(base))]  # native IDs differ and are not unique: by position
    else:
        result["aligned_by"] = None  # different number of spectra and no unique IDs: nothing to pair
        return result
    content, order_only = [], []
    for i, j in pairs:
        if base[i] == cand[j]:
            continue
        # records: [native ID, MS level, RT, peak count, ordered digest, canonical digest]
        (order_only if base[i][:4] == cand[j][:4] and base[i][5] == cand[j][5] else content).append((i, j))
    result["n_content_differs"], result["n_peak_order_differs"] = len(content), len(order_only)
    result["first"] = []
    for i, j in (content or order_only)[:limit]:  # spectra whose content differs are the more interesting ones
        entry = {"index": i, "content_differs": bool(content)}
        if j != i:
            entry["index_cand"] = j
        for position, name in enumerate(("native_id", "ms_level", "rt", "peaks")):
            entry[name] = base[i][position] if base[i][position] == cand[j][position] \
                else [base[i][position], cand[j][position]]
        result["first"].append(entry)
    return result


def compare_fingerprints(base, cand, ignore_paths=False, max_spectra=5):
    """Verdict for two fingerprints (baseline, candidate) as returned by fingerprint_experiment().

    Returns {"verdict": "IDENTICAL" | "REORDERED" | "DIFFERENT", "components": differing digests,
    "reordering": what kind of re-ordering, "count_deltas": {count: [baseline, candidate]},
    "spectra": the first differing spectra, "notes": [...]}.
    """
    result = {"verdict": "IDENTICAL", "components": [], "reordering": None, "count_deltas": {}, "notes": []}
    if base.get("fingerprint_version") != cand.get("fingerprint_version"):
        result["verdict"] = "DIFFERENT"
        result["components"] = ["fingerprint_version"]
        result["notes"].append("fingerprints were made with different versions of the fingerprint: not comparable")
        return result
    names = [name for name in ORDER_SENSITIVE if not (ignore_paths and name == "source_paths")]
    for name in names:
        left, right = base["digests"].get(name), cand["digests"].get(name)
        if left != right:
            result["components"].append(name)
        elif left is None:
            result["notes"].append("%s: not hashed (no accessor in this pyOpenMS build)" % name)
    base_counts, cand_counts = _flatten(base["counts"]), _flatten(cand["counts"])
    for key in sorted(set(base_counts) | set(cand_counts)):
        if base_counts.get(key) != cand_counts.get(key):
            result["count_deltas"][key] = [base_counts.get(key), cand_counts.get(key)]
    if base.get("path_dependent") != cand.get("path_dependent"):
        left, right = base.get("path_dependent") or [], cand.get("path_dependent") or []
        result["notes"].append("file names/paths in the metadata differ%s: %s -> %s" % (
            " (ignored)" if ignore_paths else "", json.dumps([item for item in left if item not in right]),
            json.dumps([item for item in right if item not in left])))
    if base.get("volatile") != cand.get("volatile"):
        result["notes"].append("volatile values differ (software version, completion time, tool parameters; "
                               "never part of the verdict)")
    if not result["components"] and not result["count_deltas"]:
        return result

    left, right = base["canonical"], cand["canonical"]
    not_about_order = [name for name in result["components"] if name in ("settings", "ranges", "source_paths")]
    # counts that depend on the order only: must not turn a re-ordering into DIFFERENT
    real_count_deltas = {key: value for key, value in result["count_deltas"].items()
                         if key != "spectra_not_in_canonical_peak_order"}
    if (not real_count_deltas and not not_about_order and left["ignoring_both"] == right["ignoring_both"]
            and left["chromatograms_ignoring_order"] == right["chromatograms_ignoring_order"]):
        result["verdict"] = "REORDERED"
        kinds = []
        if left["ignoring_peak_order"] != right["ignoring_peak_order"]:
            kinds.append("order of spectra")
        if left["ignoring_spectrum_order"] != right["ignoring_spectrum_order"]:
            kinds.append("order of peaks within spectra")
        if base["digests"]["chromatograms"] != cand["digests"]["chromatograms"]:
            kinds.append("order of chromatograms")
        result["reordering"] = " and ".join(kinds) + " changed, content is the same"
    else:
        result["verdict"] = "DIFFERENT"
    result["spectra"] = _spectrum_differences(base.get("per_spectrum"), cand.get("per_spectrum"), max_spectra)
    return result


# ======================================================================================
# 5. Measuring: one reader call per process
# ======================================================================================

def _maxrss_mib(children=False):
    """Peak resident set size of this process (or of its finished children) in MiB; None where unavailable."""
    if resource is None:
        return None
    value = resource.getrusage(resource.RUSAGE_CHILDREN if children else resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is in kilobytes on Linux (and the BSDs), in bytes on macOS
    return value / (1024.0 * 1024.0) if sys.platform == "darwin" else value / 1024.0


def _cpu_times(children=False):
    """(user, system) CPU seconds used so far by this process, or by its finished children."""
    if resource is not None:
        usage = resource.getrusage(resource.RUSAGE_CHILDREN if children else resource.RUSAGE_SELF)
        return usage.ru_utime, usage.ru_stime
    if children:
        return None, None
    return time.process_time(), 0.0  # Windows: user + system of this process, not separated


def _vm_hwm_mib():
    """Linux: peak RSS of the current program image (VmHWM), an independent check of ru_maxrss."""
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) / 1024.0
    except (OSError, ValueError, IndexError):
        pass
    return None


class _StdoutToStderr:
    """Send everything written to file descriptor 1 to stderr while active.

    The readers and the C++ library print progress and warnings to stdout; the worker must keep its own stdout
    clean for the one result line. Works below Python, at the OS level.
    """

    def __enter__(self):
        sys.stdout.flush()
        self._saved = os.dup(1)
        os.dup2(2, 1)
        return self

    def __exit__(self, *exc_info):
        sys.stdout.flush()
        os.dup2(self._saved, 1)
        os.close(self._saved)
        return False


def _is_bruker_directory(path):
    return os.path.isdir(path) and any(os.path.exists(os.path.join(path, name))
                                       for name in ("analysis.tdf", "analysis.tsf"))


def _path_for_type_detection(path, tmp_dir):
    """FileHandler and the TOPP tools recognise a Bruker run by its name ending in '.d'.

    A .d directory with another name (the test fixtures) is reached through a symbolic link '<tmp>/<name>.d'.
    """
    if not _is_bruker_directory(path) or path.rstrip("/\\").lower().endswith(".d"):
        return path
    link = os.path.join(tmp_dir, os.path.basename(path.rstrip("/\\")) + ".d")
    if not os.path.lexists(link):
        try:
            os.symlink(os.path.abspath(path), link, target_is_directory=True)
        except (OSError, NotImplementedError) as error:
            raise CaseSkipped("'%s' does not end in '.d' and no symbolic link could be made: %s" % (path, error))
    return link


def _set_fields(target, fields, what):
    """Assign the fields of a case to a pyOpenMS options struct (enum fields take the member name)."""
    for name, value in fields.items():
        if not hasattr(target, name):
            raise CaseSkipped("%s has no field '%s' in this pyOpenMS build" % (what, name))
        current = getattr(target, name)
        if isinstance(current, enum.Enum):
            members = type(current).__members__
            if value not in members:
                raise CaseError("%s.%s: '%s' is not one of %s" % (what, name, value, ", ".join(members)))
            value = members[value]
        try:
            setattr(target, name, value)
        except (TypeError, ValueError):
            raise CaseError("%s.%s: cannot assign %r to a field of type %s"
                            % (what, name, value, type(current).__name__))


def _public_fields(target):
    """All data fields of an options struct with their current values: the settings in effect, defaults included."""
    fields = {}
    for name in dir(target):
        if name.startswith("_"):
            continue
        value = getattr(target, name)
        if not callable(value):
            fields[name] = _jsonable(value)
    return fields


def _peak_file_options(options, settings):
    for key, value in settings.items():
        if key not in PEAK_FILE_OPTIONS:
            raise CaseError("unknown peak_file_options key '%s' (known: %s)" % (key, ", ".join(PEAK_FILE_OPTIONS)))
        setter, kind = PEAK_FILE_OPTIONS[key]
        if not hasattr(options, setter):
            raise CaseSkipped("PeakFileOptions has no %s() in this pyOpenMS build (option '%s')" % (setter, key))
        if kind == "range":
            value = oms.DRange1(float(value[0]), float(value[1]))
        elif kind == "ints":
            value = [int(level) for level in value]
        getattr(options, setter)(value)
    return options


def _prepare(case, path, job):
    """Everything before the measured call. Returns (load, read_output, effective settings, path used).

    load() is the call that is timed; it returns the MSExperiment, or None for a TOPP tool whose output file
    read_output() loads afterwards (not timed).
    """
    reader, tmp_dir = case["reader"], job["tmp_dir"]
    if reader == "bruker":
        if not hasattr(oms, "BrukerTimsFile"):
            raise CaseSkipped("this pyOpenMS has no BrukerTimsFile (built without WITH_OPENTIMS)")
        config = oms.BrukerTimsFile.Config()
        _set_fields(config, case.get("config", {}), "BrukerTimsFile.Config")
        bruker = oms.BrukerTimsFile()
        return (lambda: bruker.load(path, config)), None, _public_fields(config), path

    if reader == "thermo":
        if not hasattr(oms, "ThermoRawFile"):
            raise CaseSkipped("this pyOpenMS has no ThermoRawFile (built without WITH_THERMO_RAW)")
        thermo = oms.ThermoRawFile()
        options = thermo.getOptions()
        _set_fields(options, case.get("options", {}), "ThermoRawFileOptions")
        thermo.setOptions(options)

        def load_thermo():
            exp = oms.MSExperiment()
            thermo.load(path, exp)
            return exp
        return load_thermo, None, _public_fields(options), path

    if reader in ("filehandler", "mzml"):
        handler = oms.FileHandler() if reader == "filehandler" else oms.MzMLFile()
        handler.setOptions(_peak_file_options(handler.getOptions(), case.get("peak_file_options", {})))
        used = _path_for_type_detection(path, tmp_dir) if reader == "filehandler" else path

        def load_file():
            exp = oms.MSExperiment()
            if reader == "filehandler":
                handler.loadExperiment(used, exp)
            else:
                handler.load(used, exp)
            return exp
        return load_file, None, dict(case.get("peak_file_options", {})), used.replace(tmp_dir, "{tmp}")

    if reader == "topp":
        if not job.get("openms_bin"):
            raise CaseSkipped("TOPP cases need --openms-bin DIR")
        used = _path_for_type_detection(path, tmp_dir)
        out = os.path.join(tmp_dir, "out" + case.get("out_suffix", ".mzML"))
        try:
            command = [part.format(bin=job["openms_bin"], input=used, out=out) for part in case["command"]]
        except (KeyError, IndexError) as error:
            raise CaseError("command template: unknown placeholder %s (known: {bin}, {input}, {out})" % error)

        def run_tool():
            try:
                finished = subprocess.run(command, stdout=sys.stderr, stderr=sys.stderr)
            except OSError as error:
                raise CaseError("cannot start %s: %s" % (command[0], error))
            if finished.returncode != 0:
                raise CaseError("%s exited with code %d" % (os.path.basename(command[0]), finished.returncode))

        def read_output():
            exp = oms.MSExperiment()
            oms.MzMLFile().load(out, exp)
            os.remove(out)
            return exp
        shown = [part.replace(tmp_dir, "{tmp}") for part in command]  # the same from run to run
        return run_tool, read_output, {"command": shown}, used.replace(tmp_dir, "{tmp}")

    raise CaseError("unknown reader '%s'" % reader)


def _count_peaks(exp):
    """Peaks in all spectra (MSExperiment.getSize() would add the chromatogram data points)."""
    get_spectrum = exp.spectrum_view if hasattr(exp, "spectrum_view") else exp.__getitem__
    return sum(get_spectrum(index).size() for index in range(exp.getNrSpectra()))


def _measure(job):
    """Run one case once in this process; returns the result dict that the worker prints."""
    case, path = job["case"], job["input_path"]
    aliases = {job["tmp_dir"]: "{tmp}", os.path.realpath(job["tmp_dir"]): "{tmp}"}
    load, read_output, effective, used_path = _prepare(case, path, job)
    child = case["reader"] == "topp"  # the work happens in a child process: take its resource usage

    rss_baseline = None if child else _maxrss_mib()
    cpu_before = _cpu_times(children=child)
    started = time.perf_counter()
    exp = load()
    wall = time.perf_counter() - started
    cpu_after = _cpu_times(children=child)
    rss_peak = _maxrss_mib(children=child)  # before fingerprinting, which needs memory of its own
    vm_hwm = None if child else _vm_hwm_mib()

    result = {"status": "OK", "wall_s": wall, "rss_baseline_mib": rss_baseline, "rss_peak_load_mib": rss_peak,
              "vm_hwm_load_mib": vm_hwm, "effective_settings": effective, "loaded_from": used_path}
    if cpu_before[0] is not None:
        result["cpu_user_s"] = cpu_after[0] - cpu_before[0]
        result["cpu_sys_s"] = cpu_after[1] - cpu_before[1]
        result["cpu_s"] = result["cpu_user_s"] + result["cpu_sys_s"]
    if exp is None and job["fingerprint"]:
        exp = read_output()
    if exp is not None:
        result["n_spectra"] = exp.getNrSpectra()
        result["n_peaks"] = _count_peaks(exp)
        result["n_chromatograms"] = exp.getNrChromatograms()
    if job["fingerprint"]:
        coverage = detect_coverage(job["pinned_coverage"]) if job.get("pinned_coverage") else None
        started = time.perf_counter()
        fingerprint = fingerprint_experiment(exp, per_spectrum=job["per_spectrum"], path_aliases=aliases,
                                             coverage=coverage)
        result["fingerprint_s"] = time.perf_counter() - started
        with open(job["fingerprint_file"], "w") as handle:
            json.dump(fingerprint, handle, separators=(",", ":"), allow_nan=False)
    return result


def _probe(job):
    """What the parent needs to know about the pyOpenMS that the workers use."""
    version = oms.VersionInfo
    libraries = []
    try:  # Linux: the shared libraries that are actually loaded, to tell which build was measured
        seen = set()
        with open("/proc/self/maps") as handle:
            for line in handle:
                name = line.rstrip("\n").split(None, 5)[-1]
                if name.startswith("/") and name not in seen and \
                        re.search(r"OpenMS|pyopenms|opentims|thermo|OpenSwath", os.path.basename(name), re.I):
                    seen.add(name)
                    status = os.stat(name)
                    libraries.append({"path": name, "size": status.st_size, "mtime": status.st_mtime})
    except OSError:
        pass
    coverage = detect_coverage(job.get("pinned_coverage"))
    return {
        "status": "OK",
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "pyopenms": {"version": oms.__version__, "file": oms.__file__},
        "openms": {"version": version.getVersion(), "revision": version.getRevision(),
                   "branch": version.getBranch(), "build_time": version.getTime()},
        "readers": {"BrukerTimsFile": hasattr(oms, "BrukerTimsFile"), "ThermoRawFile": hasattr(oms, "ThermoRawFile")},
        "libraries": sorted(libraries, key=lambda library: library["path"]),
        "fingerprint_coverage": coverage["report"],
        "rss_after_import_mib": _maxrss_mib(),
        "byteorder": sys.byteorder,
    }


def worker_main(job_file):
    """Hidden entry point `_worker JOB.json`: one measurement, one result line on stdout."""
    with open(job_file) as handle:
        job = json.load(handle)
    with _StdoutToStderr():
        try:
            _import_runtime_deps()
            result = _probe(job) if job.get("probe") else _measure(job)
        except CaseSkipped as skipped:
            result = {"status": "SKIPPED", "reason": str(skipped)}
        except (CaseError, SystemExit) as error:
            result = {"status": "ERROR", "reason": str(error)}
        except BaseException as error:  # C++ exceptions arrive as Python exceptions; report all of them
            traceback.print_exc()
            result = {"status": "ERROR", "reason": "%s: %s" % (type(error).__name__, error)}
    sys.stdout.write(RESULT_MARKER + json.dumps(result) + "\n")
    sys.stdout.flush()
    return 0


# ======================================================================================
# 6. run
# ======================================================================================

_running = set()  # worker processes that are alive, so that an interrupt can stop them
_running_lock = threading.Lock()
_stopping = False  # set by an interrupt: no worker is started any more


def _stop_process(process):
    try:
        if os.name == "posix":
            os.killpg(process.pid, 9)  # the worker leads its own process group (TOPP children included)
        else:
            process.kill()
    except OSError:
        pass


def _stop_all_workers():
    """After an interrupt: the workers have sessions of their own and do not see a Ctrl-C themselves."""
    global _stopping
    with _running_lock:
        _stopping = True
        for process in list(_running):
            _stop_process(process)


def _run_worker(job, job_dir, python, timeout):
    """Start one worker process for `job`; returns (result dict, tail of its stderr)."""
    os.makedirs(job_dir, exist_ok=True)
    job_file = os.path.join(job_dir, "job.json")
    log_file = os.path.join(job_dir, "worker.log")
    with open(job_file, "w") as handle:
        json.dump(job, handle)
    command = [python, os.path.abspath(__file__), "_worker", job_file]
    with open(log_file, "wb") as log:
        with _running_lock:
            if _stopping:
                return {"status": "ERROR", "reason": "interrupted"}, []
            # a session of its own: a timeout or an interrupt can then stop the worker together with the TOPP
            # tool it may have started
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log,
                                       start_new_session=(os.name == "posix"))
            _running.add(process)
        try:
            output, _ = process.communicate(timeout=timeout)
            timed_out = False
        except subprocess.TimeoutExpired:
            _stop_process(process)
            output, _ = process.communicate()
            timed_out = True
        finally:
            with _running_lock:
                _running.discard(process)
    with open(log_file, "rb") as log:
        log.seek(max(0, os.path.getsize(log_file) - 6000))
        tail = log.read().decode("utf-8", "replace").splitlines()[-25:]
    result = None
    for line in output.decode("utf-8", "replace").splitlines():
        if line.startswith(RESULT_MARKER):
            result = json.loads(line[len(RESULT_MARKER):])
    if timed_out:
        result = {"status": "ERROR", "reason": "timeout after %g s" % timeout}
    elif result is None:
        code = process.returncode
        how = "was killed by signal %d" % -code if code < 0 else "exited with code %d" % code
        result = {"status": "ERROR", "reason": "worker %s without a result" % how}
    return result, tail


def _summary(values):
    values = [value for value in values if value is not None]
    if not values:
        return None
    return {"values": values, "min": min(values), "median": statistics.median(values), "max": max(values)}


def _run_case(index, case, context):
    """All runs of one case. Returns the report entry; a fingerprint stays in its file until the report is written."""
    args = context["args"]
    entry = {"label": case["label"], "reader": case["reader"], "input": case["input"], "status": "OK"}
    entry.update(case_settings(case))
    if case.get("comment"):
        entry["comment"] = case["comment"]
    path = context["inputs"].get(case["input"])
    if path is None:
        return dict(entry, status="SKIPPED", reason="input '%s' not given" % case["input"])
    if case["reader"] == "topp" and not args.openms_bin:
        return dict(entry, status="SKIPPED", reason="TOPP cases need --openms-bin DIR")
    runs = []
    for run in range(1 if args.parity_only else args.repeat):
        job_dir = os.path.join(context["spool"], "case%04d_run%d" % (index, run))
        tmp_dir = os.path.join(job_dir, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)
        fingerprint_file = os.path.join(context["spool"], "fingerprint%04d.json" % index)
        job = {"case": case, "input_path": path, "tmp_dir": tmp_dir, "openms_bin": args.openms_bin,
               "fingerprint": run == 0, "per_spectrum": not args.no_per_spectrum,
               "fingerprint_file": fingerprint_file, "pinned_coverage": context["pinned_coverage"]}
        result, tail = _run_worker(job, job_dir, args.python, args.timeout)
        if args.log_dir:
            name = re.sub(r"[^A-Za-z0-9_.@+-]", "_", case["label"]) + ".run%d.log" % run
            shutil.copyfile(os.path.join(job_dir, "worker.log"), os.path.join(args.log_dir, name))
        shutil.rmtree(job_dir, ignore_errors=True)
        status = result.pop("status")
        if status != "OK":
            entry.update(status=status, reason=result.get("reason", ""))
            if status == "ERROR":
                entry["log_tail"] = tail
            break
        if run == 0:
            entry["effective_settings"] = result.pop("effective_settings")
            entry["loaded_from"] = result.pop("loaded_from")
            entry["_fingerprint_file"] = fingerprint_file
        else:
            result.pop("effective_settings", None)
            result.pop("loaded_from", None)
        runs.append(result)
    if entry["status"] != "OK":
        return entry
    entry["runs"] = runs
    for key in ("wall_s", "cpu_s", "cpu_user_s", "cpu_sys_s", "rss_peak_load_mib", "rss_baseline_mib"):
        entry[key] = _summary([run.get(key) for run in runs])
    for key in ("n_spectra", "n_peaks", "n_chromatograms"):
        seen = sorted({run[key] for run in runs if key in run})
        entry[key] = seen[0] if len(seen) == 1 else seen
        if len(seen) > 1:  # the reader returned different amounts of data from run to run
            entry["unstable"] = True
    return entry


def _available_memory_gb():
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / (1024.0 * 1024.0)
    except (OSError, ValueError, IndexError):
        pass
    return None


class _MemoryBudget:
    """With --jobs > 1: let a worker start only while its estimated memory still fits.

    A whole-file Bruker load takes several GB, a frame window or a Thermo file a few hundred MB, so a fixed
    number of workers would either waste the machine on the small cases or exhaust it on the large ones.
    """

    def __init__(self, total_gb):
        self.total = self.free = total_gb
        self.changed = threading.Condition()
        self.stopped = False

    def acquire(self, amount):
        amount = min(amount, self.total)  # a case larger than the whole budget runs alone
        with self.changed:
            while self.free < amount and not self.stopped:
                self.changed.wait(0.5)
            self.free -= amount
        return amount

    def release(self, amount):
        with self.changed:
            self.free += amount
            self.changed.notify_all()


def _memory_estimate_gb(case, path, per_whole_file_gb):
    """Rough peak memory of a worker: the full reservation for a whole Bruker run, 1 GB for everything else."""
    config = case.get("config", {})
    restricted = any(key in config for key in ("frame_id_max", "rt_max_sec"))
    return per_whole_file_gb if path and os.path.isdir(path) and not restricted else 1.0


def _describe_input(path):
    info = {"path": path, "real_path": os.path.realpath(path), "kind": "directory" if os.path.isdir(path) else "file"}
    if os.path.isdir(path):
        sizes = [os.path.getsize(os.path.join(root, name)) for root, _, names in os.walk(path) for name in names]
        info["size_bytes"], info["files"] = sum(sizes), len(sizes)
    else:
        info["size_bytes"] = os.path.getsize(path)
    return info


def _git(directory, *arguments):
    try:
        finished = subprocess.run(["git", "--no-optional-locks", "-C", directory] + list(arguments),
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return finished.stdout if finished.returncode == 0 else None


def _source_revision(pyopenms_file):
    """Best effort: git revision of the source tree behind a pyopenms that is used from its build tree.

    VersionInfo carries a revision only for builds configured from a git checkout with GIT_TRACKING; this finds
    the source directory through the CMakeCache.txt above the pyopenms package and asks git (read-only).
    """
    directory = os.path.dirname(os.path.abspath(pyopenms_file))
    for _ in range(6):
        cache = os.path.join(directory, "CMakeCache.txt")
        if os.path.isfile(cache):
            with open(cache, errors="replace") as handle:
                match = re.search(r"^CMAKE_HOME_DIRECTORY:INTERNAL=(.*)$", handle.read(), re.M)
            if not match:
                return None
            source = match.group(1).strip()
            head = _git(source, "rev-parse", "HEAD")
            if head is None:
                return {"source_dir": source}
            changes = _git(source, "diff", "HEAD") or b""
            return {"source_dir": source, "revision": head.decode().strip(), "uncommitted_changes": bool(changes),
                    # identifies the state of the working tree when it has uncommitted changes
                    "diff_sha256": hashlib.sha256(changes).hexdigest() if changes else None}
        directory = os.path.dirname(directory)
    return None


def _write_report(path, header, entries):
    """Write the report; the fingerprints are spliced in from their files one at a time.

    The parent never holds the per-spectrum records of all cases in memory: a child process starts life with
    its parent's peak RSS as its own ru_maxrss, so a parent that grows would spoil the memory figures of the
    cases that run later.
    """
    temporary = path + ".tmp"
    with open(temporary, "w") as out:
        out.write('{"tool": %s, "report_version": %d,\n"header": %s,\n"cases": [\n' % (
            json.dumps(TOOL), REPORT_VERSION, json.dumps(_jsonable(header), indent=1, sort_keys=True, allow_nan=False)))
        for position, entry in enumerate(entries):
            entry = dict(entry)
            fingerprint_file = entry.pop("_fingerprint_file", None)
            text = json.dumps(_jsonable(entry), sort_keys=True, allow_nan=False)
            if fingerprint_file and os.path.exists(fingerprint_file):
                with open(fingerprint_file) as handle:
                    text = text[:-1] + ', "fingerprint": ' + handle.read().strip() + "}"
            out.write(text + (",\n" if position + 1 < len(entries) else "\n"))
        out.write("]}\n")
    os.replace(temporary, path)


def _progress_line(done, total, entry):
    text = "[%*d/%d] %-44s %-8s" % (len(str(total)), done, total, entry["label"], entry["status"])
    if entry["status"] != "OK":
        return text + entry.get("reason", "")
    text += "%8.2f s wall" % entry["wall_s"]["median"]
    if entry.get("cpu_s"):
        text += " %8.2f s cpu" % entry["cpu_s"]["median"]
    if entry.get("rss_peak_load_mib"):
        text += " %7.0f MiB" % entry["rss_peak_load_mib"]["median"]
    if isinstance(entry.get("n_spectra"), int):
        text += "  %d spectra, %d peaks" % (entry["n_spectra"], entry["n_peaks"])
    if entry.get("unstable"):
        text += "  UNSTABLE: counts differ between runs"
    return text


def cmd_run(args):
    inputs = {}
    for name, path in (("bruker_dda", args.bruker_dda), ("bruker_dia", args.bruker_dia), ("thermo", args.thermo)):
        if path:
            inputs[name] = path
    for item in args.input or []:
        name, separator, path = item.partition("=")
        if not separator or not name or not path:
            raise SystemExit("--input expects NAME=PATH, got '%s'" % item)
        inputs[name] = path
    for name, path in inputs.items():
        if not os.path.exists(path):
            raise SystemExit("input %s: '%s' does not exist" % (name, path))
        inputs[name] = os.path.abspath(path)
    if args.openms_bin:
        if not os.path.isdir(args.openms_bin):
            raise SystemExit("--openms-bin: '%s' is not a directory" % args.openms_bin)
        args.openms_bin = os.path.abspath(args.openms_bin)
    if args.repeat < 1 or args.jobs < 1:
        raise SystemExit("--repeat and --jobs must be at least 1")
    cases = select_cases(load_cases(args.cases), args.filter, args.exclude)
    if not cases:
        raise SystemExit("no case selected")
    pinned = None
    if args.coverage_from:
        with open(args.coverage_from) as handle:
            pinned = json.load(handle)["header"]["fingerprint_coverage"]["hashed"]
    if args.log_dir:
        os.makedirs(args.log_dir, exist_ok=True)

    jobs = args.jobs
    available = _available_memory_gb()
    # the memory the workers may use together: what is available now, or (unknown) room for all of them
    budget = _MemoryBudget(max(args.mem_per_job_gb, available or jobs * args.mem_per_job_gb))

    started = time.time()
    spool = tempfile.mkdtemp(prefix=TOOL + "_", dir=args.tmp_dir)
    try:
        probe, tail = _run_worker({"probe": True, "pinned_coverage": pinned}, os.path.join(spool, "probe"),
                                  args.python, 600)
        if probe.get("status") != "OK":
            print("\n".join(tail), file=sys.stderr)
            raise SystemExit("cannot use pyopenms with %s: %s" % (args.python, probe.get("reason")))
        missing = probe["fingerprint_coverage"]["missing"]
        if missing:
            print("fingerprint: not available in this pyOpenMS and therefore not hashed: %s"
                  % "; ".join("%s: %s" % (name, ", ".join(items)) for name, items in sorted(missing.items())),
                  file=sys.stderr)
        reliable = jobs == 1
        with open(os.path.abspath(__file__), "rb") as handle:
            script_sha256 = hashlib.sha256(handle.read()).hexdigest()
        header = {
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "cpus_usable": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "memory_available_gb": available,
            "python": probe["python"],
            "python_executable": args.python,
            "numpy": probe["numpy"],
            "pyopenms": probe["pyopenms"],
            "openms": probe["openms"],
            "openms_source": _source_revision(probe["pyopenms"]["file"]),
            "libraries": probe["libraries"],
            "readers": probe["readers"],
            "byteorder": probe["byteorder"],
            "rss_after_import_mib": probe["rss_after_import_mib"],
            "environment": {name: os.environ.get(name) for name in ENVIRONMENT_VARIABLES},
            "command_line": sys.argv,
            "script_sha256": script_sha256,
            "inputs": {name: _describe_input(path) for name, path in sorted(inputs.items())},
            "openms_bin": args.openms_bin,
            "settings": {"repeat": 1 if args.parity_only else args.repeat, "jobs": jobs,
                         "parity_only": args.parity_only, "per_spectrum": not args.no_per_spectrum,
                         "timeout_s": args.timeout,
                         "timings_reliable": reliable,
                         "timings_note": None if reliable else "%d workers ran at the same time: wall and CPU "
                                         "times are not comparable with a sequential run" % jobs},
            "fingerprint_version": FINGERPRINT_VERSION,
            "fingerprint_coverage": probe["fingerprint_coverage"],
            "complete": False,
        }
        context = {"args": args, "inputs": inputs, "spool": spool, "pinned_coverage": pinned}
        entries = [None] * len(cases)
        done = [0]
        lock = threading.Lock()

        def work(index):
            case = cases[index]
            reserved = budget.acquire(_memory_estimate_gb(case, inputs.get(case["input"]), args.mem_per_job_gb))
            try:
                entry = _run_case(index, case, context)
            finally:
                budget.release(reserved)
            with lock:
                entries[index] = entry
                done[0] += 1
                print(_progress_line(done[0], len(cases), entry), flush=True)
                if entry["status"] == "ERROR" and entry.get("log_tail"):
                    print("    " + "\n    ".join(entry["log_tail"][-8:]), flush=True)

        pending = list(range(len(cases)))
        signals = []  # Ctrl-C or kill: noted here, acted on in the loop below, however often it comes

        def drain():  # one of --jobs threads; each starts worker processes one after the other
            while not budget.stopped:
                with lock:
                    if not pending:
                        return
                    index = pending.pop(0)
                work(index)

        threads = [threading.Thread(target=drain, daemon=True) for _ in range(jobs)]
        previous = {number: signal.signal(number, lambda number, frame: signals.append(number))
                    for number in (signal.SIGINT, signal.SIGTERM)}
        try:
            for thread in threads:
                thread.start()
            while any(thread.is_alive() for thread in threads) and not signals:
                time.sleep(0.2)
            interrupted = bool(signals)
            if interrupted:
                print("interrupted: stopping the workers, writing what is finished", file=sys.stderr)
                budget.stopped = True
                _stop_all_workers()
                for thread in threads:
                    thread.join(10)
        finally:
            for number, handler in previous.items():
                signal.signal(number, handler)
        # a case that was cut short by the interrupt is not part of the report
        finished = [entry for entry in entries if entry is not None
                    and not (interrupted and entry["status"] == "ERROR")]
        header["complete"] = len(finished) == len(cases) and not interrupted
        header["elapsed_s"] = round(time.time() - started, 1)
        _write_report(args.out, header, finished)
    finally:
        shutil.rmtree(spool, ignore_errors=True)

    by_status = {}
    for entry in finished:
        by_status.setdefault(entry["status"], []).append(entry)
    print("%s: %d cases in %.0f s: %s" % (args.out, len(finished), header["elapsed_s"], ", ".join(
        "%d %s" % (len(group), status) for status, group in sorted(by_status.items()))))
    reasons = {}
    for entry in by_status.get("SKIPPED", []):
        reasons.setdefault(entry["reason"], []).append(entry["label"])
    for reason, labels in sorted(reasons.items()):
        print("  skipped (%s): %d cases, e.g. %s" % (reason, len(labels), labels[0]))
    for entry in by_status.get("ERROR", []):
        print("  ERROR %s: %s" % (entry["label"], entry["reason"]))
    if any(entry.get("unstable") for entry in finished):
        print("  WARNING: spectrum or peak counts differ between runs of the same case (see 'unstable')")
    return 1 if by_status.get("ERROR") or interrupted or any(entry.get("unstable") for entry in finished) else 0


# ======================================================================================
# 7. compare
# ======================================================================================

def _load_report(path):
    try:
        with open(path) as handle:
            report = json.load(handle)
    except (OSError, ValueError) as error:
        raise SystemExit("%s: %s" % (path, error))
    if report.get("tool") != TOOL or "cases" not in report:
        raise SystemExit("%s: not a report of %s" % (path, TOOL))
    return report


def _table(rows):
    widths = [max(len(row[column]) for row in rows) for column in range(len(rows[0]))]
    lines = []
    for row in rows:
        cells = [row[0].ljust(widths[0])] + [cell.rjust(width) for cell, width in zip(row[1:], widths[1:])]
        lines.append("  ".join(cells).rstrip())
    return "\n".join(lines)


def _number(summary, key, pattern="%.2f"):
    return pattern % summary[key] if summary and summary.get(key) is not None else "-"


def _format_spectrum(entry):
    def both(value):
        if isinstance(value, list):  # differs: show every digit
            return "%s -> %s" % tuple(repr(item) if isinstance(item, float) else str(item) for item in value)
        return "%.6f" % value if isinstance(value, float) else str(value)
    where = "#%d" % entry["index"] + (" (candidate #%d)" % entry["index_cand"] if "index_cand" in entry else "")
    text = "%s '%s' MS%s RT %s" % (where, both(entry["native_id"]), both(entry["ms_level"]), both(entry["rt"]))
    if "peaks" in entry:
        text += " peaks %s" % both(entry["peaks"])
    if "content_differs" in entry:
        text += ": content differs" if entry["content_differs"] else ": same peaks in another order"
    return text


def _compare_case(base, cand, args):
    """Parity entry for one pair of report cases."""
    statuses = (base["status"], cand["status"])
    if "ERROR" in statuses:
        reason = "; ".join("%s: %s" % (side, case.get("reason")) for side, case in (("baseline", base), ("candidate", cand))
                           if case["status"] == "ERROR")
        return {"verdict": "ERROR", "detail": reason}
    if statuses == ("SKIPPED", "SKIPPED"):
        return {"verdict": "SKIPPED", "detail": base.get("reason")}
    if "SKIPPED" in statuses:
        side, case = ("baseline", base) if base["status"] == "SKIPPED" else ("candidate", cand)
        return {"verdict": "MISSING", "detail": "skipped in the %s: %s" % (side, case.get("reason"))}
    if "fingerprint" not in base or "fingerprint" not in cand:
        return {"verdict": "ERROR", "detail": "no fingerprint in the report"}
    result = compare_fingerprints(base["fingerprint"], cand["fingerprint"], ignore_paths=args.ignore_paths,
                                  max_spectra=args.max_spectra)
    counts = base["fingerprint"]["counts"]
    result["size"] = "%d spectra, %d peaks, %d chromatograms" % (counts["spectra"], counts["peaks"],
                                                              counts["chromatograms"])
    for side, case in (("baseline", base), ("candidate", cand)):
        if case.get("unstable"):
            result["unstable"] = True
            result["notes"].append("NOT REPRODUCIBLE: in the %s report the runs of this case returned different "
                                   "numbers of spectra or peaks (n_spectra %s, n_peaks %s)"
                                   % (side, case.get("n_spectra"), case.get("n_peaks")))
    if base["reader"] == cand["reader"] and base.get("effective_settings") != cand.get("effective_settings"):
        keys = sorted(set(base.get("effective_settings") or {}) | set(cand.get("effective_settings") or {}))
        changed = ["%s: %r -> %r" % (key, (base.get("effective_settings") or {}).get(key, "<absent>"),
                                     (cand.get("effective_settings") or {}).get(key, "<absent>"))
                   for key in keys
                   if (base.get("effective_settings") or {}).get(key) != (cand.get("effective_settings") or {}).get(key)]
        result["notes"].append("reader settings in effect differ%s: %s" % (
            "" if args.pair else " (a default changed, or the case was edited)", "; ".join(changed)))
    return result


def _print_parity(label, result, out):
    verdict = result["verdict"]
    line = "%-52s %-10s" % (label, verdict)
    if verdict == "IDENTICAL":
        line += result["size"]  # an IDENTICAL verdict about nothing should be visible as such
    elif verdict == "REORDERED":
        line += result["reordering"]
    elif verdict == "DIFFERENT":
        line += "differs in: " + ", ".join(result["components"] or ["counts"])
    elif result.get("detail"):
        line += result["detail"]
    print(line, file=out)
    if verdict in ("DIFFERENT", "REORDERED"):
        for key, (left, right) in sorted(result["count_deltas"].items()):
            delta = " (%+d)" % (right - left) if isinstance(left, int) and isinstance(right, int) else ""
            print("    count %s: %s -> %s%s" % (key, left, right, delta), file=out)
        spectra = result.get("spectra") or {}
        if not spectra.get("available"):
            print("    (no per-spectrum records in the reports: cannot name the spectra)", file=out)
        elif spectra["aligned_by"] is None:
            print("    spectra cannot be paired: %d vs %d spectra and native IDs are not unique"
                  % (spectra["n_base"], spectra["n_cand"]), file=out)
        else:
            print("    spectra paired by %s: %d differ in content, %d only in the order of their peaks"
                  % (spectra["aligned_by"], spectra["n_content_differs"], spectra["n_peak_order_differs"]), file=out)
            for entry in spectra["first"]:
                print("      " + _format_spectrum(entry), file=out)
            for key, title in (("only_in_base", "only in the baseline"), ("only_in_cand", "only in the candidate"),
                               ("out_of_order", "at another position")):
                count = spectra.get("n_" + key, 0)
                if count:
                    print("    %d spectra %s, first:" % (count, title), file=out)
                    for entry in spectra[key]:
                        print("      " + _format_spectrum(entry), file=out)
    for note in result.get("notes", []):
        print("    note: " + note, file=out)


def cmd_compare(args):
    base_report = _load_report(args.baseline)
    cand_report = _load_report(args.candidate) if args.candidate else base_report
    out = sys.stdout
    base_header, cand_header = base_report["header"], cand_report["header"]
    if base_header.get("fingerprint_version") != cand_header.get("fingerprint_version"):
        raise SystemExit("the reports use different fingerprint versions (%s, %s): not comparable"
                         % (base_header.get("fingerprint_version"), cand_header.get("fingerprint_version")))
    warnings = []
    if base_header["fingerprint_coverage"]["hashed"] != cand_header["fingerprint_coverage"]["hashed"]:
        changes = []
        left, right = base_header["fingerprint_coverage"]["hashed"], cand_header["fingerprint_coverage"]["hashed"]
        for name in sorted(set(left) | set(right)):
            difference = sorted(set(left.get(name, [])) ^ set(right.get(name, [])))
            if difference:
                changes.append("%s: %s" % (name, ", ".join(difference)))
        warnings.append("the two pyOpenMS builds expose different accessors, so the metadata digests cover "
                        "different fields and will differ (%s). Re-run the candidate with --coverage-from %s."
                        % ("; ".join(changes), args.baseline))
    for side, header in (("baseline", base_header), ("candidate", cand_header)):
        if not header.get("complete", True):
            warnings.append("the %s report is incomplete (the run was interrupted)" % side)
    for name in sorted(set(base_header.get("inputs", {})) & set(cand_header.get("inputs", {}))):
        left, right = base_header["inputs"][name], cand_header["inputs"][name]
        if left.get("size_bytes") != right.get("size_bytes"):
            warnings.append("input '%s' is not the same data in the two reports: %s (%s bytes) and %s (%s bytes)"
                            % (name, left.get("path"), left.get("size_bytes"), right.get("path"), right.get("size_bytes")))
    if args.candidate and base_header.get("libraries") and base_header.get("libraries") == cand_header.get("libraries"):
        warnings.append("both reports were made with the very same library files (same size and modification "
                        "time): if the candidate was meant to be a new build, it was not picked up")

    base_cases = {case["label"]: case for case in base_report["cases"]}
    cand_cases = {case["label"]: case for case in cand_report["cases"]}
    if args.pair:
        pairs = []
        for left, right in args.pair:
            if left not in base_cases:
                raise SystemExit("--pair: no case '%s' in %s" % (left, args.baseline))
            if right not in cand_cases:
                raise SystemExit("--pair: no case '%s' in %s" % (right, args.candidate or args.baseline))
            pairs.append(("%s  vs  %s" % (left, right), base_cases[left], cand_cases[right]))
        only_base = only_cand = []
    else:
        if not args.candidate:
            raise SystemExit("compare needs two reports, or one report and --pair LABEL_A LABEL_B")
        labels = [case["label"] for case in base_report["cases"] if case["label"] in cand_cases]
        only_base = [case["label"] for case in base_report["cases"] if case["label"] not in cand_cases]
        only_cand = [case["label"] for case in cand_report["cases"] if case["label"] not in base_cases]
        if args.filter:
            keep = {case["label"] for case in select_cases([{"label": label} for label in
                                                            labels + only_base + only_cand], args.filter, None)}
            labels, only_base, only_cand = ([label for label in group if label in keep]
                                            for group in (labels, only_base, only_cand))
        pairs = [(label, base_cases[label], cand_cases[label]) for label in labels]

    print("baseline:  %s  (%s, %s)" % (args.baseline, base_header.get("timestamp"), _revision_text(base_header)), file=out)
    print("candidate: %s  (%s, %s)" % (args.candidate or args.baseline, cand_header.get("timestamp"),
                                       _revision_text(cand_header)), file=out)
    for warning in warnings:
        print("WARNING: " + warning, file=out)
    print("\n== output parity ==", file=out)
    results = []
    tally = {}
    for label, base, cand in pairs:
        result = _compare_case(base, cand, args)
        result["label"] = label
        results.append(result)
        tally[result["verdict"]] = tally.get(result["verdict"], 0) + 1
        _print_parity(label, result, out)
    for title, labels in (("only in the baseline", only_base), ("only in the candidate", only_cand)):
        if labels:
            print("%d cases %s: %s" % (len(labels), title, ", ".join(labels)), file=out)

    print("\n== time and memory (baseline -> candidate) ==", file=out)
    print("wall and CPU time of the reader call in seconds (min and median over the runs); RSS: peak resident set\n"
          "size during the load in MiB (median); speed-up: baseline wall min / candidate wall min", file=out)
    rows = [["case", "wall-min base", "cand", "speed-up", "wall-med base", "cand", "cpu-med base", "cand",
             "RSS base", "cand", "RSS delta"]]
    for label, base, cand in pairs:
        if base["status"] != "OK" or cand["status"] != "OK":
            continue
        base_wall, cand_wall = base.get("wall_s"), cand.get("wall_s")
        speedup = "-"
        if base_wall and cand_wall and cand_wall["min"] > 0:
            speedup = "%.2fx" % (base_wall["min"] / cand_wall["min"])
        base_rss, cand_rss = base.get("rss_peak_load_mib"), cand.get("rss_peak_load_mib")
        delta = "%+.0f" % (cand_rss["median"] - base_rss["median"]) if base_rss and cand_rss else "-"
        rows.append([label, _number(base_wall, "min"), _number(cand_wall, "min"), speedup,
                     _number(base_wall, "median"), _number(cand_wall, "median"),
                     _number(base.get("cpu_s"), "median"), _number(cand.get("cpu_s"), "median"),
                     _number(base_rss, "median", "%.0f"), _number(cand_rss, "median", "%.0f"), delta])
    if len(rows) > 1:
        print(_table(rows), file=out)
    for side, header in (("baseline", base_header), ("candidate", cand_header)):
        settings = header.get("settings", {})
        if not settings.get("timings_reliable", True):
            print("note: %s timings are unreliable: %s" % (side, settings.get("timings_note")), file=out)
        elif settings.get("repeat") == 1:
            print("note: %s timings come from a single run per case" % side, file=out)

    accepted = {"IDENTICAL", "SKIPPED"}  # SKIPPED: skipped in both reports, nothing to compare
    if args.allow_reorder:
        accepted.add("REORDERED")
    if args.allow_missing:
        accepted.add("MISSING")
    failed = [result["label"] for result in results if result["verdict"] not in accepted or result.get("unstable")]
    missing = [] if args.allow_missing or args.pair else only_base + only_cand
    compared = sum(count for verdict, count in tally.items() if verdict in ("IDENTICAL", "REORDERED", "DIFFERENT"))
    print("\n%d cases compared: %s" % (compared, ", ".join("%d %s" % (count, verdict)
                                                           for verdict, count in sorted(tally.items())) or "none"), file=out)
    exit_code = 0
    if failed or missing or compared == 0:
        exit_code = 1
        if compared == 0:
            print("FAIL: nothing was compared", file=out)
        if failed:
            print("FAIL: not %s: %s" % ("/".join(sorted(accepted - {"SKIPPED", "MISSING"})), ", ".join(failed)), file=out)
        if missing:
            print("FAIL: %d cases are in one report only (--allow-missing accepts that)" % len(missing), file=out)
    else:
        print("PASS", file=out)
    if args.json:
        with open(args.json, "w") as handle:
            json.dump({"tool": TOOL, "baseline": args.baseline, "candidate": args.candidate or args.baseline,
                       "warnings": warnings, "cases": results, "only_in_baseline": only_base,
                       "only_in_candidate": only_cand, "tally": tally, "exit_code": exit_code},
                      handle, indent=1, sort_keys=True)
    return exit_code


def _revision_text(header):
    source = header.get("openms_source") or {}
    if source.get("revision"):
        return "source %s%s" % (source["revision"][:10], " + uncommitted changes" if source.get("uncommitted_changes") else "")
    return "OpenMS %s" % (header.get("openms") or {}).get("version")


# ======================================================================================
# 8. diff: two inputs in one process, compared with tolerances
# ======================================================================================

def _parse_spec(text):
    """'reader:path[:json-options]' -> case dict. The options are the reader's settings dict."""
    reader, separator, rest = text.partition(":")
    if not separator or reader not in ("bruker", "thermo", "filehandler", "mzml"):
        raise SystemExit("'%s': expected bruker|thermo|filehandler|mzml:PATH[:JSON-OPTIONS]" % text)
    path, options = rest, {}
    position = rest.find(":{")
    if position >= 0:
        path = rest[:position]
        try:
            options = json.loads(rest[position + 1:])
        except ValueError as error:
            raise SystemExit("'%s': options are not valid JSON: %s" % (text, error))
    if not os.path.exists(path):
        raise SystemExit("'%s': '%s' does not exist" % (text, path))
    key = {"bruker": "config", "thermo": "options"}.get(reader, "peak_file_options")
    return {"label": text, "reader": reader, "input": "diff", key: options}, os.path.abspath(path)


def _load_spec(text, tmp_dir):
    case, path = _parse_spec(text)
    with _StdoutToStderr():
        try:
            load, _, _, _ = _prepare(case, path, {"tmp_dir": tmp_dir})
            return load()
        except (CaseSkipped, CaseError) as error:
            raise SystemExit("%s: %s" % (text, error))


def _match_spectra(left, right, how, rt_tolerance):
    """Pairs of spectrum indices (a, b) plus the unmatched indices of each side."""
    n_left, n_right = left.getNrSpectra(), right.getNrSpectra()
    if how == "index":
        common = min(n_left, n_right)
        return [(i, i) for i in range(common)], list(range(common, n_left)), list(range(common, n_right))
    if how == "native_id":
        ids = {}
        for index in range(n_right):
            native_id = right.spectrum_view(index).getNativeID()
            if native_id in ids:
                raise SystemExit("native ID '%s' occurs more than once in B: use --match index or rt" % native_id)
            ids[native_id] = index
        pairs, only_left, used = [], [], set()
        for index in range(n_left):
            other = ids.get(left.spectrum_view(index).getNativeID())
            if other is None or other in used:
                only_left.append(index)
            else:
                pairs.append((index, other))
                used.add(other)
        return pairs, only_left, [index for index in range(n_right) if index not in used]
    # by retention time: within one MS level, nearest neighbour inside the tolerance, each spectrum used once
    by_level = {}
    for index in range(n_right):
        spectrum = right.spectrum_view(index)
        by_level.setdefault(spectrum.getMSLevel(), []).append((spectrum.getRT(), index))
    for candidates in by_level.values():
        candidates.sort()
    pairs, only_left, used = [], [], set()
    for index in range(n_left):
        spectrum = left.spectrum_view(index)
        candidates = by_level.get(spectrum.getMSLevel(), [])
        rt = spectrum.getRT()
        best = None
        for other_rt, other in candidates:
            if other in used or abs(other_rt - rt) > rt_tolerance:
                continue
            if best is None or abs(other_rt - rt) < abs(best[0] - rt):
                best = (other_rt, other)
        if best is None:
            only_left.append(index)
        else:
            pairs.append((index, best[1]))
            used.add(best[1])
    return pairs, only_left, [index for index in range(n_right) if index not in used]


def _levels(exp):
    levels = {}
    for index in range(exp.getNrSpectra()):
        level = exp.spectrum_view(index).getMSLevel()
        levels[level] = levels.get(level, 0) + 1
    return levels


def cmd_diff(args):
    _import_runtime_deps()
    tmp_dir = tempfile.mkdtemp(prefix=TOOL + "_diff_")
    try:
        left, right = _load_spec(args.a, tmp_dir), _load_spec(args.b, tmp_dir)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    out = sys.stdout
    summary = {"a": args.a, "b": args.b, "tolerances": {"mz_ppm": args.mz_ppm, "intensity_rel": args.intensity_rel,
                                                         "rt_abs": args.rt_abs}, "match": args.match}
    for name, text, exp in (("A", args.a, left), ("B", args.b, right)):
        levels = _levels(exp)
        summary[name.lower() + "_spectra_by_ms_level"] = {str(level): count for level, count in sorted(levels.items())}
        print("%s: %s\n   %d spectra (%s), %d peaks, %d chromatograms" % (
            name, text, exp.getNrSpectra(), ", ".join("MS%d: %d" % item for item in sorted(levels.items())),
            _count_peaks(exp), exp.getNrChromatograms()), file=out)

    pairs, only_left, only_right = _match_spectra(left, right, args.match, args.rt_abs)
    worst = {"rt": (0.0, None), "mz": (0.0, None), "intensity": (0.0, None), "precursor_mz": (0.0, None)}
    level_mismatch, count_mismatch, over_mz, over_intensity, over_rt = [], [], [], [], []
    precursor_mismatch = []
    peaks_compared = 0
    for a, b in pairs:
        first, second = left.spectrum_view(a), right.spectrum_view(b)
        name = first.getNativeID() or "#%d" % a
        if first.getMSLevel() != second.getMSLevel():
            level_mismatch.append((name, first.getMSLevel(), second.getMSLevel()))
        deviation = abs(first.getRT() - second.getRT())
        if deviation > worst["rt"][0]:
            worst["rt"] = (deviation, name)
        if deviation > args.rt_abs:
            over_rt.append(name)
        first_precursors, second_precursors = first.getPrecursors(), second.getPrecursors()
        if len(first_precursors) != len(second_precursors):
            precursor_mismatch.append((name, "%d vs %d precursors" % (len(first_precursors), len(second_precursors))))
        for one, other in zip(first_precursors, second_precursors):
            reference = max(abs(one.getMZ()), abs(other.getMZ()))
            deviation = abs(one.getMZ() - other.getMZ()) / reference * 1e6 if reference > 0 else 0.0
            if deviation > worst["precursor_mz"][0]:
                worst["precursor_mz"] = (deviation, name)
            if deviation > args.mz_ppm:
                precursor_mismatch.append((name, "precursor m/z %.6f vs %.6f" % (one.getMZ(), other.getMZ())))
            elif one.getCharge() != other.getCharge():
                precursor_mismatch.append((name, "precursor charge %d vs %d" % (one.getCharge(), other.getCharge())))
        if first.size() != second.size():
            count_mismatch.append((name, first.size(), second.size()))
            continue
        if first.size() == 0:
            continue
        first_mz, first_intensity = first.get_peaks()
        second_mz, second_intensity = second.get_peaks()
        peaks_compared += len(first_mz)
        with np.errstate(divide="ignore", invalid="ignore"):
            reference = np.maximum(np.abs(first_mz), np.abs(second_mz))
            mz_deviation = np.where(reference > 0, np.abs(first_mz - second_mz) / reference * 1e6, 0.0).max()
            first_intensity = first_intensity.astype(np.float64)
            second_intensity = second_intensity.astype(np.float64)
            reference = np.maximum(np.abs(first_intensity), np.abs(second_intensity))
            intensity_deviation = np.where(reference > 0, np.abs(first_intensity - second_intensity) / reference, 0.0).max()
        if mz_deviation > worst["mz"][0]:
            worst["mz"] = (float(mz_deviation), name)
        if intensity_deviation > worst["intensity"][0]:
            worst["intensity"] = (float(intensity_deviation), name)
        if mz_deviation > args.mz_ppm:
            over_mz.append((float(mz_deviation), name))
        if intensity_deviation > args.intensity_rel:
            over_intensity.append((float(intensity_deviation), name))

    def name_of(exp, index):
        spectrum = exp.spectrum_view(index)
        return "%s (MS%d, RT %.3f, %d peaks)" % (spectrum.getNativeID() or "#%d" % index, spectrum.getMSLevel(),
                                                spectrum.getRT(), spectrum.size())

    show = args.max_spectra
    print("\nmatched by %s: %d pairs, %d only in A, %d only in B" % (args.match, len(pairs), len(only_left),
                                                                    len(only_right)), file=out)
    for title, exp, indices in (("only in A", left, only_left), ("only in B", right, only_right)):
        for index in indices[:show]:
            print("   %s: %s" % (title, name_of(exp, index)), file=out)
    print("MS level differs in %d pairs" % len(level_mismatch), file=out)
    for name, one, other in level_mismatch[:show]:
        print("   %s: MS%d vs MS%d" % (name, one, other), file=out)
    print("RT: max |difference| %.3g s (tolerance %g s), %d pairs outside%s" % (
        worst["rt"][0], args.rt_abs, len(over_rt), "; worst: %s" % worst["rt"][1] if worst["rt"][1] else ""), file=out)
    print("peak count: %d pairs equal, %d differ" % (len(pairs) - len(count_mismatch), len(count_mismatch)), file=out)
    for name, one, other in sorted(count_mismatch, key=lambda item: -abs(item[1] - item[2]))[:show]:
        print("   %s: %d vs %d peaks (%+d)" % (name, one, other, other - one), file=out)
    print("pairs with equal peak count: %d peaks compared in stored order" % peaks_compared, file=out)
    print("   m/z: max deviation %.4g ppm (tolerance %g ppm), %d spectra outside%s" % (
        worst["mz"][0], args.mz_ppm, len(over_mz), "; worst: %s" % worst["mz"][1] if worst["mz"][1] else ""), file=out)
    for deviation, name in sorted(over_mz, reverse=True)[:show]:
        print("      %s: %.4g ppm" % (name, deviation), file=out)
    print("   intensity: max relative deviation %.3g (tolerance %g), %d spectra outside%s" % (
        worst["intensity"][0], args.intensity_rel, len(over_intensity),
        "; worst: %s" % worst["intensity"][1] if worst["intensity"][1] else ""), file=out)
    for deviation, name in sorted(over_intensity, reverse=True)[:show]:
        print("      %s: %.3g" % (name, deviation), file=out)
    print("precursors: max m/z deviation %.4g ppm, %d spectra with another precursor count, m/z or charge%s" % (
        worst["precursor_mz"][0], len(precursor_mismatch),
        "" if args.check_precursors else " (reported only; --check-precursors makes it count)"), file=out)
    for name, what in precursor_mismatch[:show]:
        print("   %s: %s" % (name, what), file=out)

    problems = []
    if only_left or only_right:
        problems.append("%d unmatched spectra" % (len(only_left) + len(only_right)))
    if level_mismatch:
        problems.append("%d MS level mismatches" % len(level_mismatch))
    if over_rt:
        problems.append("%d RT differences above tolerance" % len(over_rt))
    if count_mismatch:
        problems.append("%d peak count differences" % len(count_mismatch))
    if over_mz:
        problems.append("%d spectra above the m/z tolerance" % len(over_mz))
    if over_intensity:
        problems.append("%d spectra above the intensity tolerance" % len(over_intensity))
    if precursor_mismatch and args.check_precursors:
        problems.append("%d precursor mismatches" % len(precursor_mismatch))
    if not pairs:
        problems.append("no spectra matched")
    print("\n" + ("OUTSIDE TOLERANCES: " + "; ".join(problems) if problems else "WITHIN TOLERANCES"), file=out)
    if args.json:
        summary.update({
            "matched": len(pairs), "only_in_a": len(only_left), "only_in_b": len(only_right),
            "ms_level_mismatches": len(level_mismatch), "peak_count_mismatches": len(count_mismatch),
            "peaks_compared": peaks_compared, "max_rt_difference_s": worst["rt"][0],
            "max_mz_deviation_ppm": worst["mz"][0], "max_intensity_deviation_rel": worst["intensity"][0],
            "max_precursor_mz_deviation_ppm": worst["precursor_mz"][0],
            "spectra_above_mz_tolerance": len(over_mz), "spectra_above_intensity_tolerance": len(over_intensity),
            "precursor_mismatches": len(precursor_mismatch), "problems": problems,
        })
        with open(args.json, "w") as handle:
            json.dump(summary, handle, indent=1, sort_keys=True)
    return 1 if problems else 0


# ======================================================================================
# 9. list and the command line
# ======================================================================================

def cmd_list(args):
    cases = select_cases(load_cases(args.cases), args.filter, args.exclude)
    if args.json:
        json.dump(cases, sys.stdout, indent=1)
        print()
        return 0
    rows = [["label", "reader", "input", "settings"]]
    for case in cases:
        rows.append([case["label"], case["reader"], case["input"], json.dumps(case_settings(case), sort_keys=True)])
    widths = [max(len(row[column]) for row in rows) for column in range(3)]
    for row, case in zip(rows, [None] + cases):
        print("  ".join(cell.ljust(width) for cell, width in zip(row, widths)) + "  " + row[3])
        if args.verbose and case and case.get("comment"):
            print("      " + case["comment"])
    print("%d cases" % len(cases), file=sys.stderr)
    return 0


def build_parser():
    summary, _, details = __doc__.partition("\n\n")
    parser = argparse.ArgumentParser(description=summary, epilog=details,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", metavar="{list,run,compare,diff}")

    def matrix_options(sub):
        sub.add_argument("--cases", metavar="FILE.json",
                         help="case matrix to use instead of the built-in one (a JSON list of cases; "
                              "`list --json` prints the built-in matrix in that format)")
        sub.add_argument("--filter", metavar="REGEX", help="only cases whose label matches (re.search)")
        sub.add_argument("--exclude", metavar="REGEX", help="leave out cases whose label matches")

    sub = commands.add_parser("list", help="print the case matrix",
                              description="Print the case matrix: label, reader, input, settings.")
    matrix_options(sub)
    sub.add_argument("-v", "--verbose", action="store_true", help="also say what each case exercises")
    sub.add_argument("--json", action="store_true", help="print the cases as JSON (usable with --cases)")
    sub.set_defaults(run=cmd_list)

    sub = commands.add_parser(
        "run", help="execute cases and write a JSON report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Execute cases and write a JSON report with timings, memory and output fingerprints.\n\n"
                    "Every measurement runs in a process of its own (peak RSS is a per-process high-water\n"
                    "mark). Timed is the reader call alone; the fingerprint is taken afterwards. Cases whose\n"
                    "input is not given, or that need something this pyOpenMS build does not have, are\n"
                    "recorded as SKIPPED with the reason.")
    sub.add_argument("-o", "--out", required=True, metavar="REPORT.json", help="report to write")
    sub.add_argument("--bruker-dda", metavar="PATH", help="ddaPASEF .d directory (input 'bruker_dda')")
    sub.add_argument("--bruker-dia", metavar="PATH", help="diaPASEF .d directory (input 'bruker_dia')")
    sub.add_argument("--thermo", metavar="PATH", help="Thermo .raw file (input 'thermo')")
    sub.add_argument("--input", action="append", metavar="NAME=PATH",
                     help="any other named input used by a case (repeatable)")
    sub.add_argument("--openms-bin", metavar="DIR", help="directory with the TOPP tools; without it the "
                                                         "'topp' cases are skipped")
    matrix_options(sub)
    sub.add_argument("--repeat", type=int, default=1, metavar="N",
                     help="timed runs per case, each in its own process; the fingerprint is taken in the first "
                          "(default: 1). The report has every value plus min and median")
    sub.add_argument("--parity-only", action="store_true",
                     help="one run per case, whatever --repeat says: for the parity loop, where timings do not matter")
    sub.add_argument("--jobs", type=int, default=1, metavar="N",
                     help="workers to run at the same time (default: 1). Timings are then marked as unreliable "
                          "in the report; meant for --parity-only")
    sub.add_argument("--mem-per-job-gb", type=float, default=8.0, metavar="GB",
                     help="with --jobs > 1: memory to reserve for a worker that loads a whole Bruker run (default: 8; "
                          "the default load needs 3.5 GB, MS1 aggregation more). Workers wait while the "
                          "reservations exceed the memory available at start; small cases reserve 1 GB")
    sub.add_argument("--timeout", type=float, default=3600.0, metavar="SECONDS",
                     help="per worker; a worker that takes longer is stopped and the case recorded as ERROR "
                          "(default: 3600)")
    sub.add_argument("--no-per-spectrum", action="store_true",
                     help="leave the per-spectrum records out of the report (smaller, but compare can then not "
                          "name the spectra that differ)")
    sub.add_argument("--coverage-from", metavar="REPORT.json",
                     help="hash only the metadata accessors that this earlier report hashed, so that a pyOpenMS "
                          "with more bindings stays comparable with it")
    sub.add_argument("--log-dir", metavar="DIR", help="keep everything the workers and readers print, one file per run")
    sub.add_argument("--tmp-dir", metavar="DIR", help="where to put temporary files (TOPP output can be several GB)")
    sub.add_argument("--python", default=sys.executable, metavar="EXE",
                     help="interpreter for the workers (default: the one running this script)")
    sub.set_defaults(run=cmd_run)

    sub = commands.add_parser(
        "compare", help="compare two reports: parity verdict per case, timing and memory table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Compare two reports (baseline, candidate).\n\n"
                    "Verdict per case label present in both:\n"
                    "  IDENTICAL  all digests and counts are equal\n"
                    "  REORDERED  same content, but the order of spectra, of peaks within spectra or of\n"
                    "             chromatograms changed (which of them is said)\n"
                    "  DIFFERENT  the differing components, count changes and first differing spectra are listed\n"
                    "  ERROR      the case failed in one of the runs\n"
                    "  SKIPPED    skipped in both runs (nothing to compare; does not fail)\n"
                    "  MISSING    skipped in one run only\n"
                    "Exit code 0 only if every compared case is IDENTICAL and no case is ERROR, MISSING or\n"
                    "present in only one report; 1 otherwise.\n\n"
                    "File names and paths found in the metadata form a component of their own (source_paths)\n"
                    "and count like any other; --ignore-paths leaves it out. Software versions, completion\n"
                    "times and tool parameters are reported when they differ but never change a verdict.")
    sub.add_argument("baseline", metavar="BASELINE.json")
    sub.add_argument("candidate", metavar="CANDIDATE.json", nargs="?",
                     help="may be omitted with --pair to compare two cases of one report")
    sub.add_argument("--allow-reorder", action="store_true", help="accept REORDERED as well")
    sub.add_argument("--allow-missing", action="store_true",
                     help="accept cases that are in one report only or were skipped in one run only")
    sub.add_argument("--ignore-paths", action="store_true",
                     help="ignore the component source_paths (file names and paths in the metadata), for reports "
                          "taken from the same data under another name or through a symbolic link")
    sub.add_argument("--pair", nargs=2, action="append", metavar=("LABEL_A", "LABEL_B"),
                     help="compare case LABEL_A of the baseline with case LABEL_B of the candidate (or of the same "
                          "report) instead of equal labels; repeatable. For code paths that must agree")
    sub.add_argument("--filter", metavar="REGEX", help="only cases whose label matches")
    sub.add_argument("--max-spectra", type=int, default=5, metavar="N",
                     help="differing spectra to name per case (default: 5)")
    sub.add_argument("--json", metavar="OUT.json", help="also write the comparison as JSON")
    sub.set_defaults(run=cmd_compare)

    sub = commands.add_parser(
        "diff", help="load two inputs and compare them spectrum by spectrum with tolerances",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Load two inputs in this process and compare them spectrum by spectrum with tolerances.\n\n"
                    "For readers whose output cannot be bit-identical, where hashes are useless. An input is\n"
                    "READER:PATH[:JSON-OPTIONS] with READER one of bruker, thermo, filehandler, mzml; the\n"
                    "options are BrukerTimsFile.Config fields, ThermoRawFileOptions fields or PeakFileOptions\n"
                    "keys (ms_levels, rt_range, ...).\n\n"
                    "In-process Thermo reader against an mzML from ThermoRawFileParser (doc/thermo_raw_mzml.md):\n"
                    "  * match the peak picking first: TRFP centroids by default, ThermoRawFile only with\n"
                    "    {\"centroid\": true} (or run TRFP with its no-peak-picking switch against the default);\n"
                    "    the same goes for charge, noise and all-detector exports\n"
                    "  * OpenMS keeps intensities as float32 whatever the mzML precision, so both sides are\n"
                    "    rounded alike; m/z is float64. The defaults (1 ppm, 1e-6 relative) leave room for\n"
                    "    one float32 rounding step\n"
                    "  * retention times are in seconds once loaded (TRFP writes minutes), equal to about 1e-9 s\n"
                    "  * native IDs agree (controllerType=0 controllerNumber=1 scan=N), so --match native_id works\n"
                    "  * OpenMS writes its own software/data-processing records and more vendor metadata\n"
                    "    (trailers, methods); metadata is not compared here\n"
                    "  * Precursor.getMZ() is the selected ion (the monoisotopic m/z if the instrument determined\n"
                    "    one), the isolation target is in the meta value 'isolation window target m/z'\n\n"
                    "Exit code 0 only if all spectra are matched, MS levels and peak counts agree and the\n"
                    "deviations are within the tolerances.")
    sub.add_argument("--a", required=True, metavar="READER:PATH[:JSON]", help="first input")
    sub.add_argument("--b", required=True, metavar="READER:PATH[:JSON]", help="second input")
    sub.add_argument("--mz-ppm", type=float, default=1.0, metavar="X", help="m/z tolerance in ppm (default: 1)")
    sub.add_argument("--intensity-rel", type=float, default=1e-6, metavar="X",
                     help="relative intensity tolerance |a-b|/max(|a|,|b|) (default: 1e-6)")
    sub.add_argument("--rt-abs", type=float, default=0.01, metavar="X",
                     help="retention time tolerance in seconds, also the window for --match rt (default: 0.01)")
    sub.add_argument("--match", choices=("native_id", "index", "rt"), default="native_id",
                     help="how to pair the spectra (default: native_id)")
    sub.add_argument("--check-precursors", action="store_true",
                     help="precursor count, m/z (within --mz-ppm) and charge must agree as well")
    sub.add_argument("--max-spectra", type=int, default=5, metavar="N", help="worst cases to name (default: 5)")
    sub.add_argument("--json", metavar="OUT.json", help="also write the summary as JSON")
    sub.set_defaults(run=cmd_diff)
    return parser


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) == 2 and argv[0] == "_worker":  # internal: one measurement in a fresh process
        return worker_main(argv[1])
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 2
    return args.run(args)


if __name__ == "__main__":
    sys.exit(main())
