import unittest
import os

import pyopenms

class TestSpectrumAccessOpenMS(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML").encode()

    def test_readfile_content(self):
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, exp)
        saccess = pyopenms.SpectrumAccessOpenMS(exp)
        spectrum = saccess.getSpectrumById(0)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        self.assertAlmostEqual(mz[0], 350.0000305)
        self.assertAlmostEqual(intensity[0], 0.0)
        self.assertAlmostEqual(mz[10], 358.075134277)
        self.assertAlmostEqual(intensity[10], 9210.931640625)

class TestSpectrumAccessOpenMSInMemory(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML").encode()

    def test_readfile_content(self):
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, exp)

        # Create in memory access
        saccess_ = pyopenms.SpectrumAccessOpenMS(exp)
        saccess = pyopenms.SpectrumAccessOpenMSInMemory(saccess_)

        spectrum = saccess.getSpectrumById(0)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        self.assertAlmostEqual(mz[0], 350.0000305)
        self.assertAlmostEqual(intensity[0], 0.0)
        self.assertAlmostEqual(mz[10], 358.075134277)
        self.assertAlmostEqual(intensity[10], 9210.931640625)

class TestSpectrumAccessSwathMap(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML").encode()

    def test_swathmap_openms(self):
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, exp)
        saccess = pyopenms.SpectrumAccessOpenMS(exp)

        swmap = pyopenms.SwathMap()
        swmap.lower = 400
        swmap.center = 412.5
        swmap.upper = 425
        swmap.ms1 = False

        data = swmap.getSpectrumPtr()
        assert data is None

        # Now we should have data in there
        swmap.setSpectrumPtr(saccess)
        data = swmap.getSpectrumPtr()
        assert data is not None

        swmap.setSpectrumPtr(saccess)

        data = swmap.getSpectrumPtr()
        assert data is not None

        spectrum = data.getSpectrumById(0)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        self.assertAlmostEqual(mz[0], 350.0000305)
        self.assertAlmostEqual(intensity[0], 0.0)
        self.assertAlmostEqual(mz[10], 358.075134277)
        self.assertAlmostEqual(intensity[10], 9210.931640625)

    def test_readfile_contentInMemory(self):
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, exp)

        # Create in memory access
        saccess_ = pyopenms.SpectrumAccessOpenMS(exp)
        saccess = pyopenms.SpectrumAccessOpenMSInMemory(saccess_)

        swmap = pyopenms.SwathMap()
        swmap.lower = 400
        swmap.center = 412.5
        swmap.upper = 425
        swmap.ms1 = False

        data = swmap.getSpectrumPtr()
        assert data is None

        swmap.setSpectrumPtr(saccess_)

        data = swmap.getSpectrumPtr()
        assert data is not None

        spectrum = data.getSpectrumById(0)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        self.assertAlmostEqual(mz[0], 350.0000305)
        self.assertAlmostEqual(intensity[0], 0.0)
        self.assertAlmostEqual(mz[10], 358.075134277)
        self.assertAlmostEqual(intensity[10], 9210.931640625)

if __name__ == '__main__':
    unittest.main()
