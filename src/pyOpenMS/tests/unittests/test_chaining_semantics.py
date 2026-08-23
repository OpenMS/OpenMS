"""Where a chained read or write lands, and where it silently does not (issue #9792).

Two kinds of accessor exist side by side, and a chain mixes them freely:

* **Container access** -- ``exp[i]``, iteration, ``at``/``front``/``back``,
  ``getSpectrum(i)``, and every getter returning a *collection*
  (``getPrecursors()``, ``getHits()``, ``getFloatDataArrays()``) -- returns
  **owned copies**. Writes through them do not reach the container.
* **Single sub-object accessors** -- ``getInstrumentSettings()``,
  ``getAcquisitionInfo()``, ``getSourceFile()`` -- still return a **live
  reference** into whatever object you are holding. Writes through them land in
  that object.

So a write lands exactly when the chain from a variable you hold to the target
crosses no container access. Reads are always safe: nanobind keeps the parent
alive for as long as the returned object lives, including when the parent is a
temporary produced mid-chain.

These tests pin that boundary so it cannot drift silently.
"""

import pytest


# --------------------------------------------------------------------------
# Reads are safe everywhere, including through a temporary
# --------------------------------------------------------------------------

def test_read_chain_through_container_is_safe():
    """exp[0] is a temporary; the settings alias into it and must stay readable."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    # the temporary from exp[0] is kept alive by the returned settings object
    assert exp[0].getInstrumentSettings().getZoomScan() in (True, False)
    assert exp[0].getRT() == pytest.approx(1.0)


# --------------------------------------------------------------------------
# Writes: land on a held object, lost through a container access
# --------------------------------------------------------------------------

def test_member_chain_on_held_object_lands():
    """No container access in the chain, so the write reaches the spectrum."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    spec.getInstrumentSettings().setZoomScan(True)
    assert spec.getInstrumentSettings().getZoomScan() is True


def test_member_chain_through_container_is_lost():
    """exp[0] is a copy, so the settings edited are the copy's and are discarded."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.getInstrumentSettings().setZoomScan(False)
    exp.addSpectrum(spec)

    exp[0].getInstrumentSettings().setZoomScan(True)
    assert exp[0].getInstrumentSettings().getZoomScan() is False, (
        "a write through a container access must not reach the container"
    )


def test_collection_getter_chain_is_lost():
    """getPrecursors() returns copies, so editing one does not touch the spectrum."""
    from pyopenms import MSSpectrum, Precursor

    spec = MSSpectrum()
    p = Precursor()
    p.setMZ(500.0)
    spec.setPrecursors([p])

    spec.getPrecursors()[0].setMZ(999.0)
    assert spec.getPrecursors()[0].getMZ() == pytest.approx(500.0)


def test_iteration_write_is_lost():
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
# The supported route: hold it, edit it, put it back
# --------------------------------------------------------------------------

def test_write_back_lands_including_nested_edits():
    """One write-back carries both the scalar and the nested sub-object edit."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    spec.getInstrumentSettings().setZoomScan(False)
    exp.addSpectrum(spec)

    s = exp[0]                                   # a copy we now hold
    s.setRT(77.0)
    s.getInstrumentSettings().setZoomScan(True)  # lands in the copy
    exp[0] = s                                   # and the copy replaces the element

    assert exp[0].getRT() == pytest.approx(77.0)
    assert exp[0].getInstrumentSettings().getZoomScan() is True


def test_write_back_over_iteration():
    """The index form is what replaces `for s in exp: s.setRT(...)`."""
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
