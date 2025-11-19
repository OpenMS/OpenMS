import unittest
import os

import pyopenms

class TestMobilogram(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))

    def testMobilogram(self):
        mobilogram = pyopenms.Mobilogram()
        p = pyopenms.MobilityPeak1D()
        p.setMobility(25.0)
        p.setIntensity(1e5)
        mobilogram.push_back(p)

        p_back, = list(mobilogram)
        assert isinstance(p_back, pyopenms.MobilityPeak1D)
        assert p_back.getMobility() == 25.0
        assert p_back.getIntensity() == 1e5

        mobilogram.updateRanges()
        assert isinstance(mobilogram.getMinMobility(), float)
        assert isinstance(mobilogram.getMaxMobility(), float)
        assert isinstance(mobilogram.getMinIntensity(), float)
        assert isinstance(mobilogram.getMaxIntensity(), float)

        assert mobilogram.getMinIntensity() == 1e5
        assert mobilogram.getMaxIntensity() == 1e5

    def testMobilogramRT(self):
        mobilogram = pyopenms.Mobilogram()
        mobilogram.setRT(123.45)
        assert mobilogram.getRT() == 123.45

    def testMobilogramDriftTimeUnit(self):
        mobilogram = pyopenms.Mobilogram()
        mobilogram.setDriftTimeUnit(pyopenms.DriftTimeUnit.MILLISECOND)
        assert mobilogram.getDriftTimeUnit() == pyopenms.DriftTimeUnit.MILLISECOND
        assert isinstance(mobilogram.getDriftTimeUnitAsString(), str)

    def testMobilogramGetSetPeaks(self):
        mobilogram = pyopenms.Mobilogram()
        # Set peaks using numpy arrays
        import numpy as np
        mobilities = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        intensities = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        mobilogram.set_peaks((mobilities, intensities))
        
        # Get peaks back
        mb_out, int_out = mobilogram.get_peaks()
        assert len(mb_out) == 3
        assert len(int_out) == 3
        assert mb_out[0] == 10.0
        assert mb_out[1] == 20.0
        assert mb_out[2] == 30.0
        assert int_out[0] == 100.0
        assert int_out[1] == 200.0
        assert int_out[2] == 150.0

    def testMobilogramSorting(self):
        mobilogram = pyopenms.Mobilogram()
        p1 = pyopenms.MobilityPeak1D()
        p1.setMobility(30.0)
        p1.setIntensity(100.0)
        mobilogram.push_back(p1)
        
        p2 = pyopenms.MobilityPeak1D()
        p2.setMobility(10.0)
        p2.setIntensity(200.0)
        mobilogram.push_back(p2)
        
        p3 = pyopenms.MobilityPeak1D()
        p3.setMobility(20.0)
        p3.setIntensity(150.0)
        mobilogram.push_back(p3)
        
        # Sort by position
        mobilogram.sortByPosition()
        assert mobilogram[0].getMobility() == 10.0
        assert mobilogram[1].getMobility() == 20.0
        assert mobilogram[2].getMobility() == 30.0
        assert mobilogram.isSorted()
        
        # Sort by intensity
        mobilogram.sortByIntensity(False)
        assert mobilogram[0].getIntensity() == 100.0
        assert mobilogram[1].getIntensity() == 150.0
        assert mobilogram[2].getIntensity() == 200.0

    def testMobilogramDataArrays(self):
        mobilogram = pyopenms.Mobilogram()
        fda = pyopenms.FloatDataArray()
        fda.setName("Signal to Noise Array")
        fda.push_back(15.0)
        
        mobilogram.setFloatDataArrays([fda])
        fdas = mobilogram.getFloatDataArrays()
        assert len(fdas) == 1
        assert fdas[0].getName() == b"Signal to Noise Array"

    def testMobilogramCalculateTIC(self):
        mobilogram = pyopenms.Mobilogram()
        p1 = pyopenms.MobilityPeak1D()
        p1.setMobility(10.0)
        p1.setIntensity(100.0)
        mobilogram.push_back(p1)
        
        p2 = pyopenms.MobilityPeak1D()
        p2.setMobility(20.0)
        p2.setIntensity(200.0)
        mobilogram.push_back(p2)
        
        tic = mobilogram.calculateTIC()
        assert tic == 300.0

if __name__ == '__main__':
    unittest.main()
