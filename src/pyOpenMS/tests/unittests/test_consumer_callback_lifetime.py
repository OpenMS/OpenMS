"""Regression tests for consumer-callback object lifetime (issue #9792, R-4.7 / V.5).

``NanobindMSDataConsumer`` bridges a Python consumer object to the C++
``IMSDataConsumer`` interface, so ``MzMLFile.transform``, ``MzXMLFile.transform`` and
``ImzMLFile.load(filename, consumer)`` can call back into Python.

The spectra and chromatograms it delivers live in the parser's batch buffer, which is
cleared and refilled as soon as the callback returns.  Passing them with
``rv_policy::reference`` therefore left any retained wrapper -- and every child alias
derived from it, including the zero-copy ndarray from ``get_peaks_struct()`` -- pointing
at storage that had since been reused::

    kept = []
    class C:
        def consumeSpectrum(self, s): kept.append(s)
    MzMLFile().transform(path, C())
    kept[0].getRT()          # read released storage

Copying into the callback is not a valid fix: ``IMSDataConsumer`` explicitly permits a
consumer to modify the spectra and ``MzMLHandler`` passes the modified value onward, so a
copy would silently discard the consumer's edits.  Neither is invalidating the wrapper
when the callback returns, because that cannot reach child views such as the ndarray
above.

The callback is now handed an object that *owns* its value: it is moved out of the
caller's storage for the duration of the call and moved back afterwards (so modification
still works), or copied back if the consumer kept a reference (so Python is left with a
valid independent snapshot).  These tests pin down both halves of that.
"""

import gc
import os
import sys

import numpy as np
import pytest

import pyopenms


def _data(name):
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, name)
    if not os.path.exists(path):
        pytest.skip("test data %s not available" % name)
    return path


def _churn():
    """Force collection and aggressively reuse freed heap blocks."""
    gc.collect()
    junk = [pyopenms.MSSpectrum() for _ in range(20000)]
    for s in junk:
        s.setRT(-1.0)
    arrays = [np.full(20000, -7.7) for _ in range(200)]
    del junk, arrays
    gc.collect()


class _Recorder:
    """Minimal consumer that satisfies the duck-typed interface."""

    def __init__(self):
        self.spectra = []
        self.chromatograms = []
        self.rts = []
        self.sizes = []
        self.settings_seen = 0

    def consumeSpectrum(self, s):
        self.rts.append(s.getRT())
        self.sizes.append(s.size())
        self.spectra.append(s)

    def consumeChromatogram(self, c):
        self.chromatograms.append(c)

    def setExpectedSize(self, n_spectra, n_chromatograms):
        pass

    def setExperimentalSettings(self, exp):
        self.settings_seen += 1


# ---------------------------------------------------------------------------
# A retained spectrum must stay valid
# ---------------------------------------------------------------------------


def test_retained_spectra_survive_the_transform():
    rec = _Recorder()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), rec)
    assert rec.spectra, "no spectra delivered"

    expected_rts, expected_sizes = list(rec.rts), list(rec.sizes)
    _churn()

    assert [s.getRT() for s in rec.spectra] == expected_rts
    assert [s.size() for s in rec.spectra] == expected_sizes


def test_retained_spectrum_owns_its_peak_data():
    """Not just the scalars -- the peak buffer must be the retained object's own."""
    rec = _Recorder()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), rec)
    with_peaks = [s for s in rec.spectra if s.size() > 0]
    assert with_peaks, "no spectrum with peaks delivered"

    spec = with_peaks[0]
    mz, intensity = spec.get_peaks()
    expected = (float(mz[0]), float(intensity[0]), len(mz))

    _churn()

    mz2, intensity2 = spec.get_peaks()
    assert (float(mz2[0]), float(intensity2[0]), len(mz2)) == expected


def test_retained_zero_copy_view_stays_valid():
    """A view taken inside the callback keep-alives onto the delivered object.  With an
    aliasing hand-off it pointed into the parser's buffer; it must now point into storage
    the Python object owns."""
    views = []

    class ViewKeeper(_Recorder):
        def consumeSpectrum(self, s):
            super().consumeSpectrum(s)
            if s.size() and not views:
                views.append(s.get_peaks_struct())

    rec = ViewKeeper()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), rec)
    assert views, "no view captured"

    view = views[0]
    expected = (float(view["mz"][0]), float(view["intensity"][0]), len(view))

    _churn()

    assert (float(view["mz"][0]), float(view["intensity"][0]), len(view)) == expected


def test_view_retained_without_the_spectrum_stays_valid():
    """Isolates the child keep-alive edge: the consumer keeps *only* the ndarray and drops
    every reference to the spectrum itself, so the view is the sole thing holding the
    delivered object alive."""
    views = []

    class ViewOnlyKeeper:
        def consumeSpectrum(self, s):
            if s.size() and not views:
                views.append(s.get_peaks_struct())   # note: `s` itself is not retained

        def consumeChromatogram(self, c):
            pass

        def setExpectedSize(self, a, b):
            pass

        def setExperimentalSettings(self, e):
            pass

    pyopenms.MzMLFile().transform(_data("test2.mzML"), ViewOnlyKeeper())
    assert views, "no view captured"

    view = views[0]
    expected = (float(view["mz"][0]), float(view["intensity"][0]), len(view))

    _churn()

    assert (float(view["mz"][0]), float(view["intensity"][0]), len(view)) == expected


def test_retained_chromatograms_survive_the_transform():
    """consumeChromatogram uses the same hand-off and needs its own coverage."""
    rec = _Recorder()
    pyopenms.MzMLFile().transform(_data("test.indexed.mzML"), rec)
    if not rec.chromatograms:
        pytest.skip("no chromatograms delivered by the test file")

    expected = [(c.getName(), c.size()) for c in rec.chromatograms]
    _churn()
    assert [(c.getName(), c.size()) for c in rec.chromatograms] == expected


def test_view_taken_in_callback_really_is_zero_copy():
    """Guards the premise of the test above: if get_peaks_struct() ever started copying,
    that test would pass for the wrong reason."""
    shared = []

    class ViewKeeper(_Recorder):
        def consumeSpectrum(self, s):
            super().consumeSpectrum(s)
            if s.size() and not shared:
                shared.append(bool(np.shares_memory(s.get_peaks_struct(), s.get_peaks_struct())))

    pyopenms.MzMLFile().transform(_data("test2.mzML"), ViewKeeper())
    assert shared and shared[0], "get_peaks_struct() is no longer a zero-copy view"


# ---------------------------------------------------------------------------
# The delivered object must still be a faithful, complete spectrum
# ---------------------------------------------------------------------------


def test_delivered_spectra_match_a_normal_load():
    """The move hand-off must not lose or alter anything relative to a plain load."""
    path = _data("test2.mzML")

    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(path, exp)

    rec = _Recorder()
    pyopenms.MzMLFile().transform(path, rec)

    assert len(rec.spectra) == exp.size()
    for streamed, loaded in zip(rec.spectra, exp):
        assert streamed.getRT() == loaded.getRT()
        assert streamed.size() == loaded.size()
        assert streamed.getMSLevel() == loaded.getMSLevel()
        if loaded.size():
            a_mz, a_int = streamed.get_peaks()
            b_mz, b_int = loaded.get_peaks()
            assert np.array_equal(a_mz, b_mz)
            assert np.array_equal(a_int, b_int)


def test_spectrum_is_mutable_inside_the_callback():
    """IMSDataConsumer permits a consumer to modify the spectra; the hand-off must not
    turn the argument read-only or detach it before the callback runs.

    Note this only shows the *Python* object carries the edit -- it would pass equally if
    the callback were handed a copy.  ``test_consumer_edits_reach_the_collected_map``
    below is the one that proves the edit reaches C++.
    """

    class Mutator(_Recorder):
        def consumeSpectrum(self, s):
            s.setRT(s.getRT() + 1000.0)
            s.setMSLevel(7)
            super().consumeSpectrum(s)

    rec = Mutator()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), rec)
    assert rec.spectra
    for s in rec.spectra:
        assert s.getMSLevel() == 7
    assert all(rt > 1000.0 for rt in rec.rts)


# ---------------------------------------------------------------------------
# Write-through: the consumer's edits must reach the C++ side
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "reader, path",
    [
        (pyopenms.MzMLFile, "test2.mzML"),
        (pyopenms.MzXMLFile, "test2.mzXML"),
    ],
    ids=["mzML", "mzXML"],
)
def test_consumer_edits_reach_the_collected_map(reader, path):
    """``transform(filename, consumer, exp)`` sets alwaysAppendData, so each spectrum is
    appended *after* the consumer has seen it (``MzMLHandler.cpp:186-190``).  That makes
    the write-through contract of IMSDataConsumer -- "it may modify the spectra" --
    observable from Python, and this is the test that would fail if the aliasing hand-off
    were ever "fixed" by handing the callback a copy instead.
    """

    class Mutator:
        def consumeSpectrum(self, s):
            s.setMSLevel(7)
            s.setRT(s.getRT() + 1000.0)

        def consumeChromatogram(self, c):
            pass

        def setExpectedSize(self, a, b):
            pass

        def setExperimentalSettings(self, e):
            pass

    reference = pyopenms.MSExperiment()
    reader().load(_data(path), reference)
    assert reference.size() > 0

    collected = pyopenms.MSExperiment()
    reader().transform(_data(path), Mutator(), collected)

    assert collected.size() == reference.size()
    for edited, original in zip(collected, reference):
        assert edited.getMSLevel() == 7, "consumer's edit did not reach the collected map"
        assert edited.getRT() == pytest.approx(original.getRT() + 1000.0)


def test_collected_map_is_otherwise_unchanged():
    """A consumer that touches nothing must yield exactly what a plain load gives."""

    class Noop:
        def consumeSpectrum(self, s):
            pass

        def consumeChromatogram(self, c):
            pass

        def setExpectedSize(self, a, b):
            pass

        def setExperimentalSettings(self, e):
            pass

    reference = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(_data("test2.mzML"), reference)

    collected = pyopenms.MSExperiment()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), Noop(), collected)

    assert collected.size() == reference.size()
    for got, want in zip(collected, reference):
        assert got.getRT() == want.getRT()
        assert got.size() == want.size()
        if want.size():
            a_mz, a_int = got.get_peaks()
            b_mz, b_int = want.get_peaks()
            assert np.array_equal(a_mz, b_mz)
            assert np.array_equal(a_int, b_int)


def test_retained_spectrum_is_independent_of_the_collected_map():
    """The consumer keeps its own reference while the map also collects the spectrum.
    Since Python was handed an owning object, the two must not alias: mutating the
    retained one afterwards must not disturb what was appended."""
    rec = _Recorder()
    collected = pyopenms.MSExperiment()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), rec, collected)
    assert rec.spectra and collected.size() == len(rec.spectra)

    before = [s.getRT() for s in collected]
    for s in rec.spectra:
        s.setRT(-12345.0)
    _churn()

    assert [s.getRT() for s in collected] == before


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


def test_exception_in_callback_propagates():
    class Boom(_Recorder):
        def consumeSpectrum(self, s):
            raise ValueError("a distinctive message")

    with pytest.raises(Exception) as excinfo:
        pyopenms.MzMLFile().transform(_data("test2.mzML"), Boom())

    # The bridge currently reports the failure as a RuntimeError naming the callback; the
    # original ValueError type does not survive the boundary (tracked separately). Pin the
    # behaviour down so a change is deliberate rather than accidental.
    assert "consumeSpectrum" in str(excinfo.value)


def test_exception_in_callback_does_not_corrupt_later_use():
    """The value is moved out of the parser's storage for the duration of the call, so the
    restore has to run while an exception is propagating too."""

    class BoomOnce(_Recorder):
        def __init__(self):
            super().__init__()
            self.calls = 0

        def consumeSpectrum(self, s):
            self.calls += 1
            super().consumeSpectrum(s)
            if self.calls == 1:
                raise ValueError("boom")

    rec = BoomOnce()
    with pytest.raises(Exception):
        pyopenms.MzMLFile().transform(_data("test2.mzML"), rec)

    _churn()
    # whatever was delivered before the raise must still be readable
    assert rec.spectra
    assert rec.spectra[0].getRT() == rec.rts[0]
    assert rec.spectra[0].size() == rec.sizes[0]


def test_consumer_missing_a_method_is_reported():
    class Incomplete:
        def setExpectedSize(self, a, b):
            pass

    with pytest.raises(Exception):
        pyopenms.MzMLFile().transform(_data("test2.mzML"), Incomplete())


# ---------------------------------------------------------------------------
# The ordinary streaming path must be unaffected
# ---------------------------------------------------------------------------


def test_consumer_that_discards_everything_still_works():
    class Counter:
        def __init__(self):
            self.n = 0

        def consumeSpectrum(self, s):
            self.n += 1

        def consumeChromatogram(self, c):
            pass

        def setExpectedSize(self, a, b):
            pass

        def setExperimentalSettings(self, e):
            pass

    c = Counter()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), c)
    assert c.n > 0


def test_experimental_settings_still_delivered():
    rec = _Recorder()
    pyopenms.MzMLFile().transform(_data("test2.mzML"), rec)
    assert rec.settings_seen >= 1


def test_mzxml_transform_uses_the_same_bridge():
    rec = _Recorder()
    pyopenms.MzXMLFile().transform(_data("test2.mzXML"), rec)
    assert rec.spectra

    expected = list(rec.rts)
    _churn()
    assert [s.getRT() for s in rec.spectra] == expected


# ---------------------------------------------------------------------------
# The copy branch: a consumer that both edits *and* retains
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "reader, path",
    [
        (pyopenms.MzMLFile, "test2.mzML"),
        (pyopenms.MzXMLFile, "test2.mzXML"),
    ],
    ids=["mzML", "mzXML"],
)
def test_consumer_that_edits_and_retains_gets_both(reader, path):
    """Both halves of ``restore()`` in a single consumer.

    Retaining the spectrum makes it non-unique, so ``restore()`` takes the *copy* branch:
    Python keeps the object it was handed and the caller's storage receives a copy of it.
    Every other write-through test uses a non-retaining consumer (the move branch), and
    every other retain test uses a non-editing consumer, so the intersection -- the branch
    that has to copy an edit rather than move it -- is otherwise never exercised.
    """

    class EditAndRetain:
        def __init__(self):
            self.kept = []

        def consumeSpectrum(self, s):
            s.setMSLevel(7)
            s.setRT(s.getRT() + 1000.0)
            self.kept.append(s)          # retained -> forces the copy branch

        def consumeChromatogram(self, c):
            pass

        def setExpectedSize(self, a, b):
            pass

        def setExperimentalSettings(self, e):
            pass

    reference = pyopenms.MSExperiment()
    reader().load(_data(path), reference)
    assert reference.size() > 0

    rec = EditAndRetain()
    collected = pyopenms.MSExperiment()
    reader().transform(_data(path), rec, collected)

    # the edit reached C++ through the copy
    assert collected.size() == reference.size()
    for edited, original in zip(collected, reference):
        assert edited.getMSLevel() == 7, "edit did not survive the copy branch"
        assert edited.getRT() == pytest.approx(original.getRT() + 1000.0)

    # and what Python kept is still valid, and independent of the map
    _churn()
    assert len(rec.kept) == reference.size()
    assert all(s.getMSLevel() == 7 for s in rec.kept)
    for s in rec.kept:
        s.setRT(-12345.0)
    assert all(s.getRT() != -12345.0 for s in collected), "retained object aliases the map"


def test_consumer_edits_to_chromatograms_reach_the_collected_map():
    """``consumeChromatogram`` uses the same hand-off as ``consumeSpectrum``, so its
    write-through needs its own coverage -- the spectrum tests cannot catch a regression
    that only affects chromatograms."""

    class ChromMutator:
        def __init__(self):
            self.seen = 0

        def consumeSpectrum(self, s):
            pass

        def consumeChromatogram(self, c):
            c.setName("EDITED")
            self.seen += 1

        def setExpectedSize(self, a, b):
            pass

        def setExperimentalSettings(self, e):
            pass

    rec = ChromMutator()
    collected = pyopenms.MSExperiment()
    pyopenms.MzMLFile().transform(_data("test.indexed.mzML"), rec, collected)
    if not rec.seen:
        pytest.skip("no chromatograms delivered by the test file")

    chroms = collected.getChromatograms()
    assert len(chroms) == rec.seen
    assert all(c.getName() == "EDITED" for c in chroms)


# ---------------------------------------------------------------------------
# imzML — the third caller of the bridge
# ---------------------------------------------------------------------------


def _imzml_path():
    # tests/unittests is a package, so the helper is not importable by bare name unless its
    # directory is on sys.path; test_ImzMLFile.py uses the same idiom. Without it these
    # tests would silently skip whenever they run before that file.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from test_data_paths import get_class_test_data_dir

    path = os.path.join(get_class_test_data_dir(), "ImzMLFile_1_Example_Continuous.imzML")
    if not os.path.isfile(path):
        pytest.skip("imzML test data not available: %s" % path)
    return path


def test_imzml_retained_spectra_survive_the_load():
    """``ImzMLFile.load(filename, consumer)`` is the third caller of the bridge and the only
    one that releases the GIL around the parse, so it is the only path where the callback's
    ``gil_scoped_acquire`` really acquires rather than being a re-entrant no-op. It also
    reaches the consumer through an intermediate C++ consumer (``ImzMLInterceptConsumer``),
    which forwards the very same object.
    """
    rec = _Recorder()
    pyopenms.ImzMLFile().load(_imzml_path(), rec)
    assert rec.spectra, "no spectra delivered"

    expected_rts, expected_sizes = list(rec.rts), list(rec.sizes)
    _churn()

    assert [s.getRT() for s in rec.spectra] == expected_rts
    assert [s.size() for s in rec.spectra] == expected_sizes


def test_imzml_retained_peak_data_and_views_stay_valid():
    """imzML peak data is decoded from the companion .ibd by the intermediate consumer, so
    the buffer itself is worth checking, not only the scalars."""
    views = []

    class ViewKeeper(_Recorder):
        def consumeSpectrum(self, s):
            super().consumeSpectrum(s)
            if s.size() and not views:
                views.append(s.get_peaks_struct())

    rec = ViewKeeper()
    pyopenms.ImzMLFile().load(_imzml_path(), rec)
    with_peaks = [s for s in rec.spectra if s.size() > 0]
    if not with_peaks:
        pytest.skip("no imzML spectrum with peaks delivered")

    spec = with_peaks[0]
    mz, intensity = spec.get_peaks()
    expected = (float(mz[0]), float(intensity[0]), len(mz))
    view_expected = (float(views[0]["mz"][0]), len(views[0])) if views else None

    _churn()

    mz2, intensity2 = spec.get_peaks()
    assert (float(mz2[0]), float(intensity2[0]), len(mz2)) == expected
    if view_expected is not None:
        assert (float(views[0]["mz"][0]), len(views[0])) == view_expected
