import unittest

import numpy as np
import pyopenms


def _build_im_experiment():
    exp = pyopenms.MSExperiment()

    for rt in (10.0, 20.0, 30.0):
        spec = pyopenms.MSSpectrum()
        spec.setRT(rt)
        spec.setMSLevel(2)

        precursor = pyopenms.Precursor()
        precursor.setMZ(412.5)
        precursor.setIsolationWindowLowerOffset(12.5)
        precursor.setIsolationWindowUpperOffset(12.5)
        spec.setPrecursors([precursor])

        mz = np.array([499.95, 500.00, 500.04, 510.00], dtype=np.float64)
        intensity = np.array([5.0, 100.0 + rt, 50.0 + rt, 1.0], dtype=np.float64)
        ion_mobility = np.array([0.80, 1.00, 1.03, 1.50], dtype=np.float32)

        fda = pyopenms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(ion_mobility)

        spec.set_peaks((mz, intensity))
        spec.setFloatDataArrays([fda])
        spec.sortByPosition()
        exp.addSpectrum(spec)

    exp.sortSpectra()
    return exp


def _make_coordinate():
    coord = pyopenms.ExtractionCoordinates()
    coord.mz = 500.0
    coord.rt_start = 5.0
    coord.rt_end = 25.0
    coord.ion_mobility = 1.0
    coord.id = b"tr1"
    return coord


class TestPeakMapExtractor(unittest.TestCase):

    def _assert_result(self, result):
        self.assertEqual(len(result), 1)

        peak_map = result[0]
        self.assertAlmostEqual(peak_map.target_mz, 500.0)
        self.assertAlmostEqual(peak_map.target_rt, 15.0)
        self.assertAlmostEqual(peak_map.target_ion_mobility, 1.0)
        self.assertAlmostEqual(peak_map.rt_start, 5.0)
        self.assertAlmostEqual(peak_map.rt_end, 25.0)

        self.assertEqual(list(peak_map.rt), [10.0, 10.0, 20.0, 20.0])
        np.testing.assert_allclose(list(peak_map.mz), [500.0, 500.04, 500.0, 500.04])
        np.testing.assert_allclose(list(peak_map.ion_mobility), [1.0, 1.03, 1.0, 1.03])
        np.testing.assert_allclose(list(peak_map.intensity), [110.0, 60.0, 120.0, 70.0])

    def test_extract_peak_maps_openms(self):
        exp = _build_im_experiment()
        input_map = pyopenms.SpectrumAccessOpenMS(exp)

        extractor = pyopenms.PeakMapExtractor()
        result = extractor.extractPeakMaps(input_map, [_make_coordinate()], 0.1, False, 0.1, b"tophat")

        self._assert_result(result)

    def test_extract_peak_maps_in_memory(self):
        exp = _build_im_experiment()
        input_map = pyopenms.SpectrumAccessOpenMSInMemory(pyopenms.SpectrumAccessOpenMS(exp))

        extractor = pyopenms.PeakMapExtractor()
        result = extractor.extractPeakMaps(input_map, [_make_coordinate()], 0.1, False, 0.1, b"tophat")

        self._assert_result(result)


if __name__ == "__main__":
    unittest.main()
