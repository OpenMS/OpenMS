import unittest
import os

import pyopenms

class TestOpenSwathDataStructures(unittest.TestCase):

    def setUp(self):
        pass

    def test_spectrum(self):
        # Interfaces.Spectrum uses camelCase (no addon file)
        spectrum = pyopenms.Interfaces.Spectrum()
        mz_exp = [1,2,3]
        int_exp = [4,5,6]

        spectrum.setMZArray(mz_exp)
        spectrum.setIntensityArray(int_exp)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        for m,e in zip(mz, mz_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

    def test_spectrum_osw(self):
        spectrum = pyopenms.OSSpectrum()
        mz_exp = [1,2,3]
        int_exp = [4,5,6]

        spectrum.set_mz_array(mz_exp)
        spectrum.set_intensity_array(int_exp)
        mz = spectrum.get_mz_array()
        intensity = spectrum.get_intensity_array()

        for m,e in zip(mz, mz_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

        # Now also check drift time, first check that there are 2 arrays and no drift time
        self.assertEqual( spectrum.get_drift_time_array(), None)
        self.assertEqual( len(spectrum.get_data_arrays()), 2)

        drift_exp = [7, 8, 9]
        da = pyopenms.OSBinaryDataArray()
        da.data = drift_exp
        da.description = b"Ion Mobility"
        arrays = spectrum.get_data_arrays()
        arrays.append(da)
        spectrum.set_data_arrays(arrays)

        self.assertEqual( len(spectrum.get_data_arrays()), 3)

        drift = spectrum.get_drift_time_array()
        for i,e in zip(drift, drift_exp):
            self.assertAlmostEqual(i,e)

        da = pyopenms.OSBinaryDataArray()
        da.data = [5, 6.88]
        da.description = b"test"
        arrays = spectrum.get_data_arrays()
        arrays.append(da)
        spectrum.set_data_arrays(arrays)
        self.assertEqual( len(spectrum.get_data_arrays()), 4)

        da = spectrum.get_data_arrays()
        data = da[3].get_data()
        self.assertEqual( len(data), 2)
        self.assertAlmostEqual(data[0], 5)
        self.assertAlmostEqual(data[1], 6.88)

    def test_spectrum_osw_memview(self):
        spectrum = pyopenms.OSSpectrum()
        mz_exp = [1,2,3]
        int_exp = [4,5,6]

        spectrum.set_mz_array(mz_exp)
        spectrum.set_intensity_array(int_exp)

        mz = spectrum.get_mz_array()

        self.assertAlmostEqual(mz[0], 1)

        mz_view = spectrum.get_mz_array_mv()

        self.assertAlmostEqual(mz_view[0], 1)

        # change a copy, nothing happens
        mz[0] = 100
        self.assertAlmostEqual(spectrum.get_mz_array()[0], 1)
        self.assertAlmostEqual(mz[0], 100)

        # change a memview, it changes the underlying data
        mz_view[0] = 200
        self.assertAlmostEqual(spectrum.get_mz_array()[0], 200)
        self.assertAlmostEqual(mz[0], 100)

        dataarr = spectrum.get_data_arrays()
        mz = dataarr[0].get_data()
        mz_view = dataarr[0].get_data_mv()
        self.assertAlmostEqual(mz[0], 200)
        self.assertAlmostEqual(mz_view[0], 200)

        # change a memview, it changes the underlying data
        mz[0] = 300
        mz_view[0] = 400
        self.assertAlmostEqual(spectrum.get_mz_array()[0], 400)
        self.assertAlmostEqual(mz[0], 300)

    def test_chromatogram_osw(self):
        chromatogram = pyopenms.OSChromatogram()
        rt_exp = [1,2,3]
        int_exp = [4,5,6]

        chromatogram.set_time_array(rt_exp)
        chromatogram.set_intensity_array(int_exp)
        time = chromatogram.get_time_array()
        intensity = chromatogram.get_intensity_array()

        for m,e in zip(time, rt_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

        # Now also check that we can add a data array, first check that there are 2 arrays and no drift time
        self.assertEqual( len(chromatogram.get_data_arrays()), 2)

        da = pyopenms.OSBinaryDataArray()
        da.data = [5, 6.88]
        da.description = b"test"
        arrays = chromatogram.get_data_arrays()
        arrays.append(da)
        chromatogram.set_data_arrays(arrays)
        self.assertEqual( len(chromatogram.get_data_arrays()), 3)

        da = chromatogram.get_data_arrays()
        data = da[2].get_data()
        self.assertEqual( len(data), 2)
        self.assertAlmostEqual(data[0], 5)
        self.assertAlmostEqual(data[1], 6.88)

    def test_chromatogram(self):
        # Interfaces.Chromatogram uses camelCase (no addon file)
        chromatogram = pyopenms.Interfaces.Chromatogram()
        rt_exp = [1,2,3]
        int_exp = [4,5,6]

        chromatogram.setTimeArray(rt_exp)
        chromatogram.setIntensityArray(int_exp)
        time = chromatogram.getTimeArray()
        intensity = chromatogram.getIntensityArray()

        for m,e in zip(time, rt_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

if __name__ == '__main__':
    unittest.main()
