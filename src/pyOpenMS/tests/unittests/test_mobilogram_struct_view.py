import gc
import numpy as np
import pyopenms as oms


def test_values_match_get_peaks():
    mob = oms.Mobilogram()

    mobility = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    intensity = np.array([10.0, 20.0, 30.0], dtype=np.float32)

    mob.set_peaks([mobility, intensity])

    arr = mob.get_peaks_struct()
    mob_copy, int_copy = mob.get_peaks()

    np.testing.assert_allclose(arr["mobility"], mob_copy)
    np.testing.assert_allclose(arr["intensity"], int_copy)


def test_zero_copy_modification():
    mob = oms.Mobilogram()

    mobility = np.array([5.0, 6.0], dtype=np.float64)
    intensity = np.array([50.0, 60.0], dtype=np.float32)

    mob.set_peaks([mobility, intensity])

    arr = mob.get_peaks_struct()

    arr["intensity"][0] = 999.0

    _, intensity_copy = mob.get_peaks()

    assert np.isclose(intensity_copy[0], 999.0)


def test_empty_mobilogram():
    mob = oms.Mobilogram()

    arr = mob.get_peaks_struct()

    assert isinstance(arr, np.ndarray)
    assert arr.size == 0
    assert arr.dtype['mobility'] == np.float64
    assert arr.dtype['intensity'] == np.float32


def test_single_peak():
    mob = oms.Mobilogram()

    mob.set_peaks([[1.5], [42.0]])

    arr = mob.get_peaks_struct()

    assert arr.size == 1
    assert np.isclose(arr["mobility"][0], 1.5)
    assert np.isclose(arr["intensity"][0], 42.0)


def test_large_dataset():
    mob = oms.Mobilogram()

    n = 100000

    mobility = np.linspace(1, 100, n)
    intensity = np.random.rand(n).astype(np.float32)

    mob.set_peaks([mobility, intensity])

    arr = mob.get_peaks_struct()

    assert arr.size == n

    mob_copy, int_copy = mob.get_peaks()

    np.testing.assert_allclose(arr["mobility"], mob_copy)
    np.testing.assert_allclose(arr["intensity"], int_copy)


def test_memory_lifetime():
    def create_view():
        mob = oms.Mobilogram()

        mobility = np.array([7.0], dtype=np.float64)
        intensity = np.array([77.0], dtype=np.float32)

        mob.set_peaks([mobility, intensity])

        return mob.get_peaks_struct()

    arr = create_view()

    gc.collect()

    assert np.isclose(arr["intensity"][0], 77.0)


def test_mutation_propagates():
    mob = oms.Mobilogram()

    mobility = np.array([3.0, 4.0], dtype=np.float64)
    intensity = np.array([30.0, 40.0], dtype=np.float32)

    mob.set_peaks([mobility, intensity])

    arr = mob.get_peaks_struct()

    arr["mobility"][1] = 123.0

    mob_copy, _ = mob.get_peaks()

    assert np.isclose(mob_copy[1], 123.0)
