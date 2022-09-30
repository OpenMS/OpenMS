import unittest
import os

import pyopenms

class TestOpenSwathDataStructures(unittest.TestCase):

    def setUp(self):
        pass

    def test_spectrum(self):
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

        spectrum.setMZArray(mz_exp)
        spectrum.setIntensityArray(int_exp)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        for m,e in zip(mz, mz_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

        # Now also check drift time, first check that there are 2 arrays and no drift time
        self.assertEqual( spectrum.getDriftTimeArray(), None)
        self.assertEqual( len(spectrum.getDataArrays()), 2)

        drift_exp = [7, 8, 9]
        da = pyopenms.OSBinaryDataArray()
        da.data = drift_exp 
        da.description = b"Ion Mobility";
        arrays = spectrum.getDataArrays()
        arrays.append(da)
        spectrum.setDataArrays(arrays)

        self.assertEqual( len(spectrum.getDataArrays()), 3)

        drift = spectrum.getDriftTimeArray()
        for i,e in zip(drift, drift_exp):
            self.assertAlmostEqual(i,e)

        da = pyopenms.OSBinaryDataArray()
        da.data = [5, 6.88] 
        da.description = b"test";
        arrays = spectrum.getDataArrays()
        arrays.append(da)
        spectrum.setDataArrays(arrays)
        self.assertEqual( len(spectrum.getDataArrays()), 4)

        da = spectrum.getDataArrays()
        data = da[3].getData()
        self.assertEqual( len(data), 2)
        self.assertAlmostEqual(data[0], 5)
        self.assertAlmostEqual(data[1], 6.88)

    def test_spectrum_osw_memview(self):
        spectrum = pyopenms.OSSpectrum()
        mz_exp = [1,2,3]
        int_exp = [4,5,6]

        spectrum.setMZArray(mz_exp)
        spectrum.setIntensityArray(int_exp)

        mz = spectrum.getMZArray()

        self.assertAlmostEqual(mz[0], 1)

        mz_view = spectrum.getMZArray_mv()

        self.assertAlmostEqual(mz_view[0], 1)

        # change a copy, nothing happens
        mz[0] = 100
        self.assertAlmostEqual(spectrum.getMZArray()[0], 1)
        self.assertAlmostEqual(mz[0], 100)

        # change a memview, it changes the underlying data
        mz_view[0] = 200
        self.assertAlmostEqual(spectrum.getMZArray()[0], 200)
        self.assertAlmostEqual(mz[0], 100)

        dataarr = spectrum.getDataArrays()
        mz = dataarr[0].getData()
        mz_view = dataarr[0].getData_mv()
        self.assertAlmostEqual(mz[0], 200)
        self.assertAlmostEqual(mz_view[0], 200)

        # change a memview, it changes the underlying data
        mz[0] = 300
        mz_view[0] = 400
        self.assertAlmostEqual(spectrum.getMZArray()[0], 400)
        self.assertAlmostEqual(mz[0], 300)

    def test_chromatogram_osw(self):
        chromatogram = pyopenms.OSChromatogram()
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

        # Now also check that we can add a data array, first check that there are 2 arrays and no drift time
        self.assertEqual( len(chromatogram.getDataArrays()), 2)

        da = pyopenms.OSBinaryDataArray()
        da.data = [5, 6.88] 
        da.description = b"test";
        arrays = chromatogram.getDataArrays()
        arrays.append(da)
        chromatogram.setDataArrays(arrays)
        self.assertEqual( len(chromatogram.getDataArrays()), 3)

        da = chromatogram.getDataArrays()
        data = da[2].getData()
        self.assertEqual( len(data), 2)
        self.assertAlmostEqual(data[0], 5)
        self.assertAlmostEqual(data[1], 6.88)

    def test_chromatogram(self):
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
