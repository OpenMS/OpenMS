"""Regression tests for set_peaks() leaving aligned metadata misaligned (issue #9792, R-6.2).

``set_peaks()`` replaces the peak table but does not touch the parallel binary data arrays
(float, string, integer).  Replacing the peaks of a container that carries such an array
therefore used to leave it silently inconsistent::

    spec.setFloatDataArrays([im_array])   # 100 entries, aligned with 100 peaks
    spec.set_peaks((mz5, int5))           # 5 peaks, array still has 100 entries

Nothing complained at the call.  The corruption surfaced much later -- as an exception from
a ``sortByPosition()``/``select()`` far from the cause, or as an mzML file written with a
mis-sized binary data array.  ``set_peaks()`` now applies a ``metadata`` policy: ``error``
(the default) refuses and leaves the container untouched, ``clear`` drops the arrays that no
longer fit, and ``keep`` is the old behaviour, opted into explicitly.

Three deliberate design points are pinned down here, because all three are easy to
"simplify" away:

* **The trigger is a size mismatch, not peak replacement.**  An array that is already the
  right length is kept, with no exception and no kwarg -- see
  ``test_same_length_replacement_keeps_arrays``.  Invalidating on every ``set_peaks()`` call
  would destroy the annotations of correct programs that edit m/z or intensity in place
  (``test_IonMobilityData.py:135`` does exactly that), and would re-do what commit 7215fc74af
  (``[FIX] do not clear data arrays when setting peaks``, bug #5876) deliberately undid.
* **The check runs against the new peak count, before the container is resized.**  Checking
  ``self.size()`` instead would compare the arrays against the count they are still aligned
  with and accept the very case this exists to catch -- see
  ``test_shrink_with_stale_array_raises`` and ``test_grow_with_stale_array_raises``.
* **The "clear" drop happens after the peaks are written, and only to the stranded arrays.**
  The incoming intensities may be a zero-copy view into the very array being dropped, so
  freeing first is a use-after-free (``test_clear_does_not_free_an_aliasing_input_array``);
  and an array that still fits is a valid annotation, not collateral damage
  (``test_clear_policy_keeps_arrays_that_still_fit``).
"""

import numpy as np
import pytest

import pyopenms


# (class, noun used in the exception message)
_CONTAINERS = [
    pytest.param(pyopenms.MSSpectrum, "spectrum", id="MSSpectrum"),
    pytest.param(pyopenms.MSChromatogram, "chromatogram", id="MSChromatogram"),
    pytest.param(pyopenms.Mobilogram, "mobilogram", id="Mobilogram"),
]

# All three spellings are in use across the existing suite, and Mobilogram's sequence overload
# is a hand-written tuple/list branch rather than the generic one the other two use.
_CALL_FORMS = [
    pytest.param(lambda c, p, i, **kw: c.set_peaks(p, i, **kw), id="two_positional"),
    pytest.param(lambda c, p, i, **kw: c.set_peaks((p, i), **kw), id="tuple"),
    pytest.param(lambda c, p, i, **kw: c.set_peaks([p, i], **kw), id="list"),
]

_ARRAY_KINDS = ["float", "string", "integer"]


def _peaks(n, descending=False):
    """n peaks with strictly monotonic positions and distinguishable intensities."""
    pos = np.arange(1.0, n + 1.0, dtype=np.float64)
    if descending:
        pos = pos[::-1].copy()
    return pos, (np.arange(1.0, n + 1.0) * 10.0).astype(np.float32)


def _attach_array(container, kind, size, name="Annotation"):
    """Attach a single named data array of `size` entries of the given family."""
    if kind == "float":
        arr = pyopenms.FloatDataArray()
        arr.set_data(np.arange(size, dtype=np.float32))
        arr.setName(name)
        container.setFloatDataArrays([arr])
    elif kind == "integer":
        arr = pyopenms.IntegerDataArray()
        arr.set_data(np.arange(size, dtype=np.int32))
        arr.setName(name)
        container.setIntegerDataArrays([arr])
    elif kind == "string":
        arr = pyopenms.StringDataArray()
        arr.set_data([str(i) for i in range(size)])
        arr.setName(name)
        container.setStringDataArrays([arr])
    else:
        raise AssertionError("unknown array kind %r" % kind)


def _arrays(container, kind):
    return {
        "float": container.getFloatDataArrays,
        "string": container.getStringDataArrays,
        "integer": container.getIntegerDataArrays,
    }[kind]()


def _make(cls, n_peaks, kind=None, n_array=None, descending=False):
    """A container with `n_peaks` peaks and, optionally, one data array of `n_array` entries."""
    container = cls()
    pos, inten = _peaks(n_peaks, descending=descending)
    container.set_peaks(pos, inten)
    if kind is not None:
        _attach_array(container, kind, n_peaks if n_array is None else n_array)
    return container


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_shrink_with_stale_array_raises(cls, noun, call, kind):
    """Replacing 100 peaks with 5 would strand a 100-entry array, so the call is refused."""
    container = _make(cls, 100, kind)
    pos, inten = _peaks(5)
    with pytest.raises(RuntimeError, match="does not match the new %s size" % noun):
        call(container, pos, inten)


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_grow_with_stale_array_raises(cls, noun, kind):
    """Growing strands the array just as shrinking does; the check is not one-sided."""
    container = _make(cls, 5, kind)
    pos, inten = _peaks(100)
    with pytest.raises(RuntimeError, match="does not match the new %s size" % noun):
        container.set_peaks(pos, inten)


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_rejected_call_leaves_container_unchanged(cls, noun, call, kind):
    """A refused set_peaks() is a strict no-op: peaks, array *values* and array name survive."""
    container = _make(cls, 10, kind)
    before_pos, before_int = container.get_peaks()
    before_data = list(_arrays(container, kind)[0].get_data())

    with pytest.raises(RuntimeError):
        call(container, *_peaks(3))

    after_pos, after_int = container.get_peaks()
    assert container.size() == 10
    assert np.array_equal(after_pos, before_pos)
    assert np.array_equal(after_int, before_int)
    arrays = _arrays(container, kind)
    assert len(arrays) == 1
    assert list(arrays[0].get_data()) == before_data
    assert arrays[0].getName() == "Annotation"


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_same_length_replacement_keeps_arrays(cls, noun, call, kind):
    """The in-place edit idiom -- same peak count, new values -- keeps annotations untouched.

    This is the #5876 guarantee and needs no ``metadata`` argument.
    """
    container = _make(cls, 4, kind)
    pos, inten = _peaks(4)
    call(container, pos + 100.0, inten)

    assert container.size() == 4
    assert np.allclose(container.get_peaks()[0], pos + 100.0)
    arrays = _arrays(container, kind)
    assert len(arrays) == 1
    assert len(arrays[0].get_data()) == 4
    assert arrays[0].getName() == "Annotation"


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("policy", ["error", "clear", "keep"])
def test_no_data_arrays_is_unaffected(cls, noun, call, policy):
    """The common case -- no annotations at all -- behaves identically under every policy."""
    container = cls()
    pos, inten = _peaks(6)
    call(container, pos, inten, metadata=policy)

    assert container.size() == 6
    assert np.allclose(container.get_peaks()[0], pos)


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_empty_named_array_is_exempt(cls, noun, kind):
    """An empty array carries no per-peak association, so it cannot be mis-associated.

    This matches the ``!arrays[i].empty()`` exemption in the core's own check; without it the
    binding would be stricter than MSChromatogram::checkDataArraySizes_().
    """
    container = _make(cls, 3, kind, n_array=0)
    container.set_peaks(*_peaks(7))

    assert container.size() == 7
    arrays = _arrays(container, kind)
    assert len(arrays) == 1
    assert len(arrays[0].get_data()) == 0
    assert arrays[0].getName() == "Annotation"


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_empty_named_array_survives_clear(cls, noun, kind):
    """An exempt array is exempt under every policy, including the destructive one."""
    container = _make(cls, 3, kind, n_array=0)
    container.set_peaks(*_peaks(7), metadata="clear")

    arrays = _arrays(container, kind)
    assert len(arrays) == 1
    assert arrays[0].getName() == "Annotation"


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
def test_clear_policy_drops_all_three_families(cls, noun, call):
    """metadata="clear" must not forget a family -- all of float, string and integer go."""
    container = _make(cls, 10)
    for kind in _ARRAY_KINDS:
        _attach_array(container, kind, 10)

    pos, inten = _peaks(2)
    call(container, pos, inten, metadata="clear")

    assert container.size() == 2
    for kind in _ARRAY_KINDS:
        assert len(_arrays(container, kind)) == 0


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_clear_policy_keeps_arrays_that_still_fit(cls, noun, call, kind):
    """"clear" drops what the new count stranded, not everything.

    The policy is documented as being about arrays the new peak count leaves mis-sized, so an
    array that still matches is not collateral damage -- it is still a valid annotation.
    """
    container = _make(cls, 4, kind)
    pos, inten = _peaks(4)
    call(container, pos + 100.0, inten, metadata="clear")

    arrays = _arrays(container, kind)
    assert len(arrays) == 1
    assert len(arrays[0].get_data()) == 4
    assert arrays[0].getName() == "Annotation"


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
def test_clear_policy_drops_only_the_stranded_array(cls, noun):
    """A mis-sized float array goes; an aligned integer array stays."""
    container = _make(cls, 6)
    _attach_array(container, "float", 6)      # aligned with the NEW count below
    _attach_array(container, "integer", 99)   # stranded

    container.set_peaks(*_peaks(6), metadata="clear")

    assert len(_arrays(container, "float")) == 1
    assert len(_arrays(container, "integer")) == 0


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("kind", _ARRAY_KINDS)
def test_keep_policy_reproduces_legacy_behaviour(cls, noun, call, kind):
    """metadata="keep" is the pre-fix behaviour, available to callers who realign themselves."""
    container = _make(cls, 10, kind)
    call(container, *_peaks(2), metadata="keep")

    assert container.size() == 2
    arrays = _arrays(container, kind)
    assert len(arrays) == 1
    assert len(arrays[0].get_data()) == 10


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
def test_clear_does_not_free_an_aliasing_input_array(cls, noun):
    """The incoming intensities may be a zero-copy view *into* the array "clear" is about to drop.

    ``FloatDataArray.get_data_view()`` is zero-copy and ``getFloatDataArrays()`` hands out the
    container's own arrays, so dropping them before the peaks are written would leave set_peaks()
    filling from released memory. The drop therefore happens after the fill.
    """
    container = _make(cls, 64, "float", n_array=64)
    aliased = _arrays(container, "float")[0].get_data_view()[:8]
    expected = np.array(aliased, dtype=np.float32)   # detached copy of what must be stored

    container.set_peaks(_peaks(8)[0], aliased, metadata="clear")

    assert container.size() == 8
    assert np.allclose(container.get_peaks()[1], expected)
    assert len(_arrays(container, "float")) == 0


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("call", _CALL_FORMS)
@pytest.mark.parametrize("policy", ["drop", "", "Error", "ignore"])
def test_invalid_policy_raises_value_error(cls, noun, call, policy):
    container = cls()
    with pytest.raises(ValueError, match="Must be 'error', 'clear' or 'keep'"):
        call(container, *_peaks(3), metadata=policy)


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
def test_invalid_policy_does_not_mutate(cls, noun):
    """The policy is parsed before any storage is touched, including under "clear"."""
    container = _make(cls, 5, "float")
    with pytest.raises(ValueError):
        container.set_peaks(*_peaks(2), metadata="drop")

    assert container.size() == 5
    assert len(_arrays(container, "float")[0].get_data()) == 5


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("policy", ["error", "clear", "keep"])
def test_input_length_mismatch_leaves_arrays_intact(cls, noun, policy):
    """Mismatched inputs are rejected before the policy runs, so "clear" cannot half-apply."""
    container = _make(cls, 5, "float")
    pos, _ = _peaks(5)
    _, short_int = _peaks(3)

    with pytest.raises(RuntimeError, match="same length"):
        container.set_peaks(pos, short_int, metadata=policy)

    assert container.size() == 5
    assert len(_arrays(container, "float")) == 1
    assert len(_arrays(container, "float")[0].get_data()) == 5


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
def test_clear_lets_a_later_sort_succeed(cls, noun):
    """The reported symptom: a mis-sized array makes a later sort throw far from the cause.

    The peaks are descending so the sort actually permutes -- an already-sorted container
    short-circuits without inspecting the arrays.
    """
    kept = _make(cls, 10, "float", descending=True)
    kept.set_peaks(*_peaks(4, descending=True), metadata="keep")
    with pytest.raises(RuntimeError):
        kept.sortByPosition()

    cleared = _make(cls, 10, "float", descending=True)
    cleared.set_peaks(*_peaks(4, descending=True), metadata="clear")
    cleared.sortByPosition()
    assert np.allclose(cleared.get_peaks()[0], _peaks(4)[0])


@pytest.mark.parametrize("cls,noun", _CONTAINERS)
@pytest.mark.parametrize("seq", [tuple, list])
@pytest.mark.parametrize("policy", ["clear", b"clear"])
def test_policy_passed_positionally_with_a_sequence_is_reported(cls, noun, seq, policy):
    """set_peaks(peaks, "clear") binds to the two-argument overload, which is registered first.

    Without an explicit guard it fails inside the numpy conversion with a message that says
    nothing about ``metadata``.  ``bytes`` is covered too, because the shared std::string caster
    accepts it, so ``set_peaks(mz, intensity, b"clear")`` is a working spelling and the
    sequence form has to diagnose it the same way.
    """
    container = cls()
    pos, inten = _peaks(3)
    with pytest.raises(ValueError, match="metadata="):
        container.set_peaks(seq((pos, inten)), policy)
