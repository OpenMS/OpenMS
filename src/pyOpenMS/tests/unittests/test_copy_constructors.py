"""
Tests that copy.copy(), copy.deepcopy(), and ClassName(obj) work for all
pyopenms classes that support copying.

Catches regressions where nanobind dispatches copy-construction to a string
constructor instead of the copy constructor due to registration order.

Uses a subprocess probe at collection time to discover classes that segfault
during default-construction or copy. Those classes are excluded from the
parametrized tests but cause test_no_segfaults_during_copy to FAIL with
a clear list of crashers.
"""

import copy
import enum
import inspect
import subprocess
import sys

import pytest

import pyopenms


def _get_copyable_classes():
    """Return public non-enum classes that define __copy__ directly (not inherited).

    Classes that merely inherit __copy__ from a base class (e.g., DefaultParamHandler,
    ProgressLogger) suffer from object slicing: copy.copy() returns the base type.
    Those classes need their own __copy__ binding before they can be tested here.
    """
    return sorted(
        name for name in dir(pyopenms)
        if not name.startswith("_")
        and inspect.isclass(getattr(pyopenms, name, None))
        and not issubclass(getattr(pyopenms, name), enum.Enum)
        and "__copy__" in getattr(pyopenms, name).__dict__
    )


_PROBE_CODE = """\
import copy, sys
import pyopenms
for name in {class_names!r}:
    cls = getattr(pyopenms, name)
    try:
        obj = cls()
    except Exception:
        print("SKIP:" + name, flush=True)
        continue
    try:
        copy.copy(obj)
    except Exception:
        print("SKIP:" + name, flush=True)
        continue
    print("OK:" + name, flush=True)
"""


def _probe_safe_classes(class_names, timeout=120):
    """Probe which classes can be default-constructed and copied without segfaulting.

    Runs operations in a subprocess. If it crashes, identifies the crasher from
    the last successfully printed class name and recurses on remaining classes.

    Returns (safe_set, crashed_list).
    """
    if not class_names:
        return set(), []

    code = _PROBE_CODE.format(class_names=class_names)
    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        return set(), list(class_names)

    safe = set()
    skipped = set()
    for line in result.stdout.strip().splitlines():
        if line.startswith("OK:"):
            safe.add(line[3:])
        elif line.startswith("SKIP:"):
            skipped.add(line[5:])

    if result.returncode == 0:
        return safe, []

    # Subprocess crashed. Identify crasher and recurse on remaining.
    reported = safe | skipped
    remaining = [n for n in class_names if n not in reported]
    crashed = []
    if remaining:
        crashed.append(remaining[0])  # first unreported class = the crasher
        remaining = remaining[1:]
    if remaining:
        more_safe, more_crashed = _probe_safe_classes(remaining, timeout)
        safe |= more_safe
        crashed.extend(more_crashed)

    return safe, crashed


_ALL_COPYABLE = _get_copyable_classes()
_SAFE_CLASSES, _CRASHED_CLASSES = _probe_safe_classes(_ALL_COPYABLE)
COPYABLE_CLASSES = sorted(_SAFE_CLASSES)


def _try_default_construct(cls):
    """Try to default-construct a class instance. Returns instance or None."""
    try:
        return cls()
    except (TypeError, RuntimeError):
        return None


def _has_own_eq(cls):
    """Check if cls defines __eq__ (not inherited from object)."""
    return any(
        "__eq__" in base.__dict__
        for base in cls.__mro__
        if base is not object
    )


def test_no_segfaults_during_copy():
    """Fail if any copyable class segfaults during default-construct + copy."""
    assert not _CRASHED_CLASSES, (
        f"{len(_CRASHED_CLASSES)} class(es) segfault during copy probe: "
        + ", ".join(_CRASHED_CLASSES)
    )


@pytest.mark.parametrize("class_name", COPYABLE_CLASSES)
def test_copy_copy(class_name):
    cls = getattr(pyopenms, class_name)
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    obj_copy = copy.copy(obj)
    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


@pytest.mark.parametrize("class_name", COPYABLE_CLASSES)
def test_deepcopy(class_name):
    cls = getattr(pyopenms, class_name)
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    obj_copy = copy.deepcopy(obj)
    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


@pytest.mark.parametrize("class_name", COPYABLE_CLASSES)
def test_copy_constructor(class_name):
    cls = getattr(pyopenms, class_name)
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    try:
        obj_copy = cls(obj)
    except TypeError:
        pytest.skip(f"{class_name} does not support copy construction")

    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


# Classes whose copy.copy() and copy.deepcopy() worked in pyOpenMS 3.5.0. Most of them
# derive from a bound base (DefaultParamHandler, ProgressLogger, CVTermList, CsvFile)
# and need their own __copy__/__deepcopy__, or they inherit the base's and return a
# copy of the base class instead of themselves (#10260).
COPYABLE_IN_3_5 = [
    "AScore", "AbsoluteQuantitation", "AbsoluteQuantitationMethodFile",
    "AccurateMassSearchEngine", "Biosaur2Algorithm", "CachedMzMLHandler", "ChannelInfo",
    "ConfidenceScoring", "Contact", "DTA2DFile", "DataFilter", "DataFilters",
    "ElutionModelFitter", "ElutionPeakDetection", "EmgGradientDescent", "EmgScoring",
    "FeatureFinderMultiplexAlgorithm", "GNPSMGFFile", "GaussFilter", "IDFilter",
    "ILPDCWrapper", "InternalCalibration", "IsotopeLabelingMDVs", "MRMAssay", "MRMDecoy",
    "MRMFeatureFilter", "MRMFeaturePickerFile", "MRMFeatureQCFile", "MRMTransitionGroupPicker",
    "MS2File", "MZTrafoModel", "MascotGenericFile", "MassTraceDetection", "MassTraces",
    "MasstraceCorrelator", "MetaboliteSpectralMatching", "MultiplexDeltaMasses",
    "MultiplexDeltaMassesGenerator", "OpenPepXLAlgorithm", "PeakIntegrator",
    "PeakPickerChromatogram", "PeakPickerHiRes", "PeakPickerIM", "PeakPickerIterative",
    "PeptideAndProteinQuant", "PeptideIndexing", "Prediction", "Protein", "Publication",
    "SavitzkyGolayFilter", "SeedListGenerator", "SimpleSearchEngineAlgorithm",
    "SiriusExportAlgorithm", "TargetedExperiment_Instrument",
    "TargetedExperiment_Interpretation", "TargetedExperiment_Modification",
    "TargetedSpectraExtractor", "TraMLProduct", "TransitionPQPFile", "TransitionTSVFile",
    "XFDRAlgorithm",
]


@pytest.mark.parametrize("class_name", COPYABLE_IN_3_5)
def test_copy_keeps_the_class(class_name):
    cls = getattr(pyopenms, class_name)
    obj = cls()
    for duplicate in (copy.copy(obj), copy.deepcopy(obj), cls(obj)):
        assert type(duplicate) is cls


# Further classes with a copy constructor that derive from DefaultParamHandler, whose
# __copy__ they inherited, so copy.copy() returned a DefaultParamHandler
DERIVED_WITH_OWN_COPY = {
    "FeatureDistance": lambda: pyopenms.FeatureDistance(1.0, False),
    "MultiplexResolverAlgorithm": lambda: pyopenms.MultiplexResolverAlgorithm(),
    "SwathMapMassCorrection": lambda: pyopenms.SwathMapMassCorrection(),
}


@pytest.mark.parametrize("class_name", sorted(DERIVED_WITH_OWN_COPY))
def test_copy_of_derived_class_keeps_the_class(class_name):
    cls = getattr(pyopenms, class_name)
    obj = DERIVED_WITH_OWN_COPY[class_name]()
    for duplicate in (copy.copy(obj), copy.deepcopy(obj), cls(obj)):
        assert type(duplicate) is cls


def _has_copy_constructor(cls):
    """True if a constructor of cls takes a single argument of type cls."""
    for signature, *_ in getattr(cls.__init__, "__nb_signature__", ()):
        params = signature[signature.index("(") + 1:signature.rindex(")")].split(", ")
        args = [p for p in params if p not in ("self", "/", "*")]
        if len(args) == 1 and args[0].split(": ")[-1].split(".")[-1] == cls.__name__:
            return True
    return False


def test_classes_with_a_copy_constructor_do_not_inherit_copy():
    # A class that inherits __copy__ from a bound base (DefaultParamHandler, ProgressLogger,
    # XMLFile, ...) gets a copy of that base from copy.copy() (#10260)
    assert _has_copy_constructor(pyopenms.MSSpectrum)  # the signature check works
    inheriting = sorted(
        name for name, cls in inspect.getmembers(pyopenms, inspect.isclass)
        if _has_copy_constructor(cls) and "__copy__" not in cls.__dict__ and hasattr(cls, "__copy__"))
    assert inheriting == []


def _peak_width_estimator():
    # PeakWidthEstimator fits a spline to the widths of picked peaks, so it needs some
    mzs = [400.0 + 5.0 * i for i in range(200)]
    spectrum = pyopenms.MSSpectrum()
    spectrum.set_peaks((mzs, [100.0] * len(mzs)))
    experiment = pyopenms.MSExperiment()
    experiment.addSpectrum(spectrum)
    boundaries = []
    for mz in mzs:
        boundary = pyopenms.PeakBoundary()
        boundary.mz_min, boundary.mz_max = mz * (1 - 1e-5), mz * (1 + 1e-5)
        boundaries.append(boundary)
    return pyopenms.PeakWidthEstimator(experiment, [boundaries])


# Classes that must refuse to be copied, each with a way to make one
UNCOPYABLE = {
    # Has no C++ copy constructor (it owns a unique_ptr<SimpleSVM>). Without its own
    # __copy__, copy.copy() returned the DefaultParamHandler part of it.
    "FeatureFindingMetabo": lambda tmp_path: pyopenms.FeatureFindingMetabo(),
    # These delete objects that their C++ copy would share, so the copy and the original
    # would both delete them. MSDataSqlConsumer(other) aborted Python in 3.5.0.
    "MSDataSqlConsumer": lambda tmp_path: pyopenms.MSDataSqlConsumer(str(tmp_path / "x.sqMass"), 1, 500, True, False, 1e-4),
    "CachedSwathFileConsumer": lambda tmp_path: pyopenms.CachedSwathFileConsumer(str(tmp_path) + "/", "cached", 0, []),
    "MzMLSwathFileConsumer": lambda tmp_path: pyopenms.MzMLSwathFileConsumer(str(tmp_path) + "/", "mzml", 0, []),
    "PeakWidthEstimator": lambda tmp_path: _peak_width_estimator(),
}


@pytest.mark.parametrize("name", sorted(UNCOPYABLE))
def test_uncopyable_class_refuses_to_copy(name, tmp_path):
    obj = UNCOPYABLE[name](tmp_path)
    with pytest.raises(TypeError, match="cannot be copied"):
        copy.copy(obj)
    with pytest.raises(TypeError, match="cannot be copied"):
        copy.deepcopy(obj)
    with pytest.raises(TypeError):
        type(obj)(obj)  # no copy constructor X(other)
