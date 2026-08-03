"""Regression tests for iterator lifetime (keep_alive on __iter__).

nanobind's ``make_iterator`` returns an object that stores *raw* C++ iterators
into the container.  The iterator holds no reference to the container itself,
and ``rv_policy::reference_internal`` only governs the elements that are
yielded -- it does not keep the container alive.  Every ``__iter__`` binding
must therefore carry ``nb::keep_alive<0, 1>()``.

Without it, the following is a use-after-free::

    it = iter(build_spectrum())   # temporary dies right here
    list(it)                      # reads freed memory

These tests build the container as a temporary whose only possible owner is
the iterator, drop every other reference, churn the heap so that a freed
buffer is very likely to be reused, and then check that iteration still
yields the original values.
"""

import gc
import sys

import pytest

import pyopenms


def _churn():
    """Force collection and aggressively reuse freed heap blocks."""
    gc.collect()
    # Allocations of assorted sizes, then dropped, to reuse a freed buffer.
    junk = [bytearray(b"\xa5" * n) for n in (16, 64, 256, 1024, 4096) for _ in range(200)]
    del junk
    gc.collect()


# --- container factories -------------------------------------------------
#
# Each entry is (name, factory, extractor).  ``factory`` returns a freshly
# built container; ``extractor`` maps a yielded element to a comparable value.


def _spectrum():
    s = pyopenms.MSSpectrum()
    for i in range(64):
        p = pyopenms.Peak1D()
        p.setMZ(100.0 + i)
        p.setIntensity(float(i))
        s.push_back(p)
    return s


def _chromatogram():
    c = pyopenms.MSChromatogram()
    for i in range(64):
        p = pyopenms.ChromatogramPeak()
        p.setRT(100.0 + i)
        p.setIntensity(float(i))
        c.push_back(p)
    return c


def _experiment():
    e = pyopenms.MSExperiment()
    for i in range(32):
        s = pyopenms.MSSpectrum()
        s.setRT(100.0 + i)
        e.addSpectrum(s)
    return e


def _mobilogram():
    m = pyopenms.Mobilogram()
    for i in range(64):
        p = pyopenms.MobilityPeak1D()
        p.setMobility(100.0 + i)
        p.setIntensity(float(i))
        m.push_back(p)
    return m


def _feature_map():
    fm = pyopenms.FeatureMap()
    for i in range(32):
        f = pyopenms.Feature()
        f.setRT(100.0 + i)
        fm.push_back(f)
    return fm


def _consensus_map():
    cm = pyopenms.ConsensusMap()
    for i in range(32):
        f = pyopenms.ConsensusFeature()
        f.setRT(100.0 + i)
        cm.push_back(f)
    return cm


def _param():
    p = pyopenms.Param()
    for i in range(32):
        p.setValue("k%02d" % i, i, "doc")
    return p


def _peptide_id_list():
    lst = pyopenms.PeptideIdentificationList()
    for i in range(32):
        pid = pyopenms.PeptideIdentification()
        pid.setRT(100.0 + i)
        lst.push_back(pid)
    return lst


def _aa_sequence():
    return pyopenms.AASequence.fromString("DFPIANGERDFPIANGERDFPIANGER")


def _isotope_distribution():
    gen = pyopenms.CoarseIsotopePatternGenerator(8)
    return pyopenms.EmpiricalFormula("C100H200N30O40S2").getIsotopeDistribution(gen)


def _na_sequence():
    return pyopenms.NASequence.fromString("AUGCAUGCAUGCAUGC")


def _peak_group():
    pg = pyopenms.PeakGroup()
    for i in range(64):
        p = pyopenms.LogMzPeak()
        p.mz = 100.0 + i
        pg.push_back(p)
    return pg


def _deconvolved_spectrum():
    ds = pyopenms.DeconvolvedSpectrum()
    for i in range(32):
        pg = pyopenms.PeakGroup()
        pg.setIndex(i)
        ds.push_back(pg)
    return ds


CONTAINERS = [
    ("MSSpectrum", _spectrum, lambda p: p.getMZ()),
    ("MSChromatogram", _chromatogram, lambda p: p.getRT()),
    ("MSExperiment", _experiment, lambda s: s.getRT()),
    ("Mobilogram", _mobilogram, lambda p: p.getMobility()),
    ("FeatureMap", _feature_map, lambda f: f.getRT()),
    ("ConsensusMap", _consensus_map, lambda f: f.getRT()),
    ("Param", _param, lambda e: e.name),
    ("PeptideIdentificationList", _peptide_id_list, lambda p: p.getRT()),
    ("AASequence", _aa_sequence, lambda r: r.getOneLetterCode()),
    ("IsotopeDistribution", _isotope_distribution, lambda p: p.getMZ()),
    ("NASequence", _na_sequence, lambda r: r.getCode()),
    ("PeakGroup", _peak_group, lambda p: p.mz),
    ("DeconvolvedSpectrum", _deconvolved_spectrum, lambda pg: pg.getIndex()),
]

IDS = [c[0] for c in CONTAINERS]


def _expected(factory, extract):
    """Reference values, read while the container is provably alive.

    This must NOT iterate a temporary: on a build without the keep-alive the
    temporary is freed before iteration, so computing the baseline that way
    would itself trigger the use-after-free under test and could compare
    garbage against garbage.
    """
    container = factory()
    values = [extract(x) for x in container]
    del container
    return values


@pytest.mark.parametrize("name,factory,extract", CONTAINERS, ids=IDS)
def test_iterator_holds_a_reference_to_its_container(name, factory, extract):
    """The keep-alive edge must be observable, not merely probable.

    The other tests in this file provoke undefined behaviour and rely on it
    manifesting.  This one checks the annotation directly: nb::keep_alive<0, 1>()
    increments the container's refcount for as long as the iterator lives, so a
    missing annotation fails here deterministically, with no UB involved.
    """
    container = factory()

    before = sys.getrefcount(container)
    it = iter(container)
    during = sys.getrefcount(container)
    assert during == before + 1, (
        "%s.__iter__ did not retain its container "
        "(refcount %d -> %d); nb::keep_alive<0, 1>() is missing" % (name, before, during)
    )

    del it
    gc.collect()
    assert sys.getrefcount(container) == before


@pytest.mark.parametrize("name,factory,extract", CONTAINERS, ids=IDS)
def test_iterator_outlives_temporary_container(name, factory, extract):
    """iter() on a temporary must keep that temporary alive."""
    expected = _expected(factory, extract)
    assert expected, "%s factory produced an empty container" % name

    # The container is a temporary: after this statement the iterator is the
    # only thing that could possibly own it.
    it = iter(factory())
    _churn()

    assert [extract(x) for x in it] == expected


@pytest.mark.parametrize("name,factory,extract", CONTAINERS, ids=IDS)
def test_iterator_outlives_deleted_container(name, factory, extract):
    """An iterator must survive `del` of the last named reference."""
    container = factory()
    expected = [extract(x) for x in container]

    it = iter(container)
    del container
    _churn()

    assert [extract(x) for x in it] == expected


@pytest.mark.parametrize("name,factory,extract", CONTAINERS, ids=IDS)
def test_partially_consumed_iterator_survives(name, factory, extract):
    """Dropping the container mid-iteration must not invalidate the rest."""
    expected = _expected(factory, extract)
    if len(expected) < 4:
        pytest.skip("%s is too small for a partial-consumption test" % name)

    container = factory()
    it = iter(container)
    seen = [extract(next(it)), extract(next(it))]

    del container
    _churn()

    seen.extend(extract(x) for x in it)
    assert seen == expected
