"""Load an mzML file into a light, numpy based run summary.

Only what the detectors need is kept: centroided peaks of the MSn spectra
(optionally subsampled), precursor information, isolation windows, activation,
filter strings, and the file-level instrument description.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

PROTON = 1.007276466812


@dataclass
class Scan:
    native_id: str
    ms_level: int
    rt: float                     # seconds
    mz: np.ndarray
    intensity: np.ndarray
    precursor_mz: float = 0.0
    charge: int = 0
    iso_lower: float = 0.0        # isolation window offsets (Th)
    iso_upper: float = 0.0
    activation: tuple = ()        # short names, e.g. ("HCD",)
    energy: float = 0.0
    filter_string: str = ""
    scan_low: float = 0.0
    scan_high: float = 0.0
    ion_mobility: bool = False
    drift_time: float = -1.0
    tic: float = 0.0

    @property
    def iso_width(self) -> float:
        return self.iso_lower + self.iso_upper

    @property
    def iso_center(self) -> float:
        # the window centre, which for DIA is not the "precursor"
        return self.precursor_mz + (self.iso_upper - self.iso_lower) / 2.0


@dataclass
class Run:
    path: str
    instrument_name: str = ""
    instrument_vendor: str = ""
    analyzers: list = field(default_factory=list)   # file-level analyzer names
    software: list = field(default_factory=list)
    ms1: list = field(default_factory=list)          # list[Scan]
    ms2: list = field(default_factory=list)
    ms3: list = field(default_factory=list)
    n_spectra: int = 0
    n_ms1_total: int = 0
    n_ms2_total: int = 0
    n_ms3_total: int = 0

    def summary(self) -> dict:
        return {
            "path": self.path,
            "instrument": self.instrument_name,
            "vendor": self.instrument_vendor,
            "analyzers": self.analyzers,
            "spectra": self.n_spectra,
            "ms1": self.n_ms1_total,
            "ms2": self.n_ms2_total,
            "ms3": self.n_ms3_total,
        }


_FILTER_CE = re.compile(r"@([a-z]+)(\d+(?:\.\d+)?)", re.I)


def _text(value) -> str:
    """pyopenms returns bytes for some strings in older builds (e.g. the Python 3.10 wheels)."""
    return value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value)


def _activation_from_filter(fs: str):
    """Thermo filter strings name every activation step: '@hcd27.00' or
    '@etd50.00@hcd25.00' (EThcD)."""
    steps = _FILTER_CE.findall(fs or "")
    return [(name.upper(), float(e)) for name, e in steps]


def load_mzml(path: str, max_ms2: int = 20000, max_ms1: int = 2000,
              seed: int = 0) -> Run:
    """Read an mzML file. MS2 spectra beyond max_ms2 are subsampled uniformly
    so that the detectors stay fast on large runs (the counts in Run keep the
    totals)."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    run = Run(path=str(path))
    ins = exp.getInstrument()
    run.instrument_name = _text(ins.getName() or ins.getModel() or "")
    run.instrument_vendor = _text(ins.getVendor() or "")
    for ma in ins.getMassAnalyzers():
        name = str(ma.getType()).split(".")[-1]
        if name != "ANALYZERNULL":
            run.analyzers.append(name)
    try:
        run.software = [s.getName() for s in ins.getSoftware()]
    except Exception:  # older pyopenms
        pass

    spectra = list(exp)
    run.n_spectra = len(spectra)
    levels = np.array([s.getMSLevel() for s in spectra])
    run.n_ms1_total = int((levels == 1).sum())
    run.n_ms2_total = int((levels == 2).sum())
    run.n_ms3_total = int((levels >= 3).sum())

    rng = np.random.default_rng(seed)

    def pick(idx, limit):
        if len(idx) <= limit:
            return idx
        return np.sort(rng.choice(idx, size=limit, replace=False))

    keep = set()
    keep.update(pick(np.where(levels == 1)[0], max_ms1).tolist())
    keep.update(pick(np.where(levels == 2)[0], max_ms2).tolist())
    keep.update(np.where(levels >= 3)[0][:max_ms2].tolist())

    for i in sorted(keep):
        s = spectra[i]
        mz, inten = s.get_peaks()
        mz = np.asarray(mz, dtype=np.float64)
        inten = np.asarray(inten, dtype=np.float64)
        fs = ""
        if s.metaValueExists("filter string"):
            fs = str(s.getMetaValue("filter string"))
        sw = s.getInstrumentSettings().getScanWindows()
        scan = Scan(
            native_id=s.getNativeID(),
            ms_level=s.getMSLevel(),
            rt=s.getRT(),
            mz=mz,
            intensity=inten,
            filter_string=fs,
            scan_low=sw[0].begin if sw else (float(mz.min()) if mz.size else 0.0),
            scan_high=sw[0].end if sw else (float(mz.max()) if mz.size else 0.0),
            ion_mobility=bool(s.containsIMData()) or s.getDriftTime() > 0,
            drift_time=s.getDriftTime(),
            tic=float(inten.sum()),
        )
        precs = s.getPrecursors()
        if precs:
            p = precs[0]
            scan.precursor_mz = p.getMZ()
            scan.charge = p.getCharge()
            scan.iso_lower = p.getIsolationWindowLowerOffset()
            scan.iso_upper = p.getIsolationWindowUpperOffset()
            scan.activation = tuple(sorted(_text(m) for m in p.getActivationMethodsAsShortString()))
            scan.energy = p.getActivationEnergy()
            if p.metaValueExists("drift time") or p.getDriftTime() > 0:
                scan.ion_mobility = True
        steps = _activation_from_filter(fs)
        if steps:
            if not scan.activation:
                scan.activation = tuple(sorted({n for n, _ in steps}))
            if scan.energy <= 0:
                scan.energy = steps[-1][1]
        if scan.ms_level == 1:
            run.ms1.append(scan)
        elif scan.ms_level == 2:
            run.ms2.append(scan)
        else:
            run.ms3.append(scan)
    return run


def write_subset_mzml(src: str, dst: str, n_top_tic: int = 3000) -> int:
    """Write the MS2 spectra with the highest TIC (and all MS1) to dst, for a
    quick pre-search. Returns the number of MS2 spectra written."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(src), exp)
    ms2 = [(s.calculateTIC(), i) for i, s in enumerate(exp) if s.getMSLevel() == 2]
    ms2.sort(reverse=True)
    chosen = {i for _, i in ms2[:n_top_tic]}
    # pyopenms has no MSExperiment.setExperimentalSettings: keep the run metadata by replacing the spectra
    keep = [s for i, s in enumerate(exp.getSpectra()) if s.getMSLevel() == 1 or i in chosen]
    exp.setSpectra(keep)
    exp.updateRanges()
    oms.MzMLFile().store(str(dst), exp)
    return len(chosen)
