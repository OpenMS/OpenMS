"""Where a chained read or write lands, and where it silently does not (issue #9792).

The rule is uniform: **every getter returns an owned copy**. A write therefore
lands only in the object you are holding, and reaches its parent only when you
put it back with the matching setter.

    s = exp[0]                                  # a copy
    settings = s.getInstrumentSettings()        # a copy of the copy
    settings.setZoomScan(True)
    s.setInstrumentSettings(settings)           # into the spectrum we hold
    exp[0] = s                                  # into the experiment

Two kinds of accessor still return live references, and both are self-evident
at the call site:

* a **fluent builder** returning itself for chaining (``ParquetFilter.eq(...).andNext()``)
* a **database lookup** returning a shared, process-lifetime entry
  (``ResidueDB``, ``ModificationsDB``)

Reads are always safe, including through a temporary produced mid-chain:
nanobind keeps the parent alive for as long as the object taken from it lives.

These tests pin the boundary so it cannot drift silently.
"""

import pytest


# --------------------------------------------------------------------------
# Reads are safe everywhere
# --------------------------------------------------------------------------

def test_read_chain_through_container_is_safe():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    assert exp[0].getInstrumentSettings().getZoomScan() is False
    assert exp[0].getRT() == pytest.approx(1.0)


# --------------------------------------------------------------------------
# Every chained write is discarded -- there is no exception to learn
# --------------------------------------------------------------------------

def test_chained_write_on_held_object_is_discarded():
    """Even with no container in the chain: getInstrumentSettings() returns a copy."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    spec.getInstrumentSettings().setZoomScan(True)
    assert spec.getInstrumentSettings().getZoomScan() is False


def test_chained_write_through_container_is_discarded():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    exp.addSpectrum(MSSpectrum())

    exp[0].getInstrumentSettings().setZoomScan(True)
    assert exp[0].getInstrumentSettings().getZoomScan() is False


def test_collection_getter_chain_is_discarded():
    from pyopenms import MSSpectrum, Precursor

    spec = MSSpectrum()
    p = Precursor()
    p.setMZ(500.0)
    spec.setPrecursors([p])

    spec.getPrecursors()[0].setMZ(999.0)
    assert spec.getPrecursors()[0].getMZ() == pytest.approx(500.0)


def test_iteration_write_is_discarded():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    for rt in (1.0, 2.0):
        s = MSSpectrum()
        s.setRT(rt)
        exp.addSpectrum(s)

    for s in exp:
        s.setRT(99.0)
    assert [exp[i].getRT() for i in range(2)] == pytest.approx([1.0, 2.0])


# --------------------------------------------------------------------------
# The supported route: hold it, edit it, put it back -- at every level
# --------------------------------------------------------------------------

def test_write_back_one_level():
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    settings = spec.getInstrumentSettings()
    settings.setZoomScan(True)
    spec.setInstrumentSettings(settings)
    assert spec.getInstrumentSettings().getZoomScan() is True


def test_write_back_two_levels():
    """A nested edit needs a write-back at each level it crossed."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    s = exp[0]
    s.setRT(77.0)
    settings = s.getInstrumentSettings()
    settings.setZoomScan(True)
    s.setInstrumentSettings(settings)
    exp[0] = s

    assert exp[0].getRT() == pytest.approx(77.0)
    assert exp[0].getInstrumentSettings().getZoomScan() is True


def test_write_back_over_iteration():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    for rt in (1.0, 2.0):
        s = MSSpectrum()
        s.setRT(rt)
        exp.addSpectrum(s)

    for i in range(len(exp)):
        s = exp[i]
        s.setRT(s.getRT() * 10.0)
        exp[i] = s

    assert [exp[i].getRT() for i in range(2)] == pytest.approx([10.0, 20.0])


def test_collection_write_back():
    from pyopenms import MSSpectrum, Precursor

    spec = MSSpectrum()
    p = Precursor()
    p.setMZ(500.0)
    spec.setPrecursors([p])

    precs = spec.getPrecursors()
    precs[0].setMZ(999.0)
    spec.setPrecursors(precs)
    assert spec.getPrecursors()[0].getMZ() == pytest.approx(999.0)


# --------------------------------------------------------------------------
# The two documented exceptions
# --------------------------------------------------------------------------

def test_fluent_builder_still_chains():
    """A builder returns itself, so chaining is the whole point and must keep working."""
    from pyopenms import ParquetFilter

    f = ParquetFilter()
    assert f.eq("ms_level", 1).andNext().eq("charge", 2) is not None


def test_struct_attributes_are_live_not_copies():
    """Plain data-holder structs expose members as attributes, which behave as
    Python attributes do: they are the object, so edits land. This is the one
    place the copy rule does not apply, and it is pinned so it cannot drift."""
    from pyopenms import AQS_featureConcentration

    fc = AQS_featureConcentration()
    fc.feature.setRT(42.0)
    assert fc.feature.getRT() == pytest.approx(42.0), (
        "attribute access must stay live; only getters return copies"
    )
