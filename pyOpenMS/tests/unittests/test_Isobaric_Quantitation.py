import unittest
import os

import pyopenms

def general_setup(self):
    self.eightplex = pyopenms.ItraqEightPlexQuantitationMethod()
    self.fourplex = pyopenms.ItraqFourPlexQuantitationMethod()
    self.tmt = pyopenms.TMTSixPlexQuantitationMethod()

class TestIsobaricQuantitationMethod(unittest.TestCase):

    def setUp(self):
        pass

    def test_init_fxn(self):
        eightplex = pyopenms.ItraqEightPlexQuantitationMethod()
        fourplex = pyopenms.ItraqFourPlexQuantitationMethod()
        tmt = pyopenms.TMTSixPlexQuantitationMethod()

        for inst in [eightplex, fourplex, tmt]:
            assert inst.getName() is not None
            assert inst.getChannelInformation() is not None
            assert inst.getNumberOfChannels() is not None
            assert inst.getIsotopeCorrectionMatrix() is not None
            assert inst.getReferenceChannel() is not None

        self.assertEqual(tmt.getIsotopeCorrectionMatrix().rows(), 6)
        self.assertEqual(tmt.getIsotopeCorrectionMatrix().cols(), 6)

        self.assertEqual(fourplex.getIsotopeCorrectionMatrix().rows(), 4)
        self.assertEqual(fourplex.getIsotopeCorrectionMatrix().cols(), 4)

        self.assertEqual(eightplex.getIsotopeCorrectionMatrix().rows(), 8)
        self.assertEqual(eightplex.getIsotopeCorrectionMatrix().cols(), 8)

class TestIsobaricChannelExtractor(unittest.TestCase):

    def setUp(self):
        general_setup(self)

    def testInit(self):
        assert pyopenms.IsobaricChannelExtractor(self.eightplex) is not None
        assert pyopenms.IsobaricChannelExtractor(self.fourplex) is not None
        assert pyopenms.IsobaricChannelExtractor(self.tmt) is not None

    def testFunction(self):

        for method in [self.eightplex, self.fourplex, self.tmt]:
            inst = pyopenms.IsobaricChannelExtractor(method)
            map1 = pyopenms.ConsensusMap()
            exp = pyopenms.MSExperiment()
            assert inst.extractChannels is not None
            # Will not work with empty map
            # inst.extractChannels(exp, map1)

class TestIsobaricNormalizer(unittest.TestCase):

    def setUp(self):
        general_setup(self)

    def testInit(self):
        assert pyopenms.IsobaricNormalizer(self.eightplex) is not None
        assert pyopenms.IsobaricNormalizer(self.fourplex) is not None
        assert pyopenms.IsobaricNormalizer(self.tmt) is not None

    def testFunction(self):

        for method in [self.eightplex, self.fourplex, self.tmt]:
            inst = pyopenms.IsobaricNormalizer(method)
            map1 = pyopenms.ConsensusMap()
            assert inst.normalize is not None
            inst.normalize(map1)

class TestIsobaricIsotopeCorrector(unittest.TestCase):

    def setUp(self):
        general_setup(self)

    def testInit(self):
        assert pyopenms.IsobaricIsotopeCorrector(self.eightplex) is not None
        assert pyopenms.IsobaricIsotopeCorrector(self.fourplex) is not None
        assert pyopenms.IsobaricIsotopeCorrector(self.tmt) is not None

    def testFunction(self):

        for method in [self.eightplex, self.fourplex, self.tmt]:
            inst = pyopenms.IsobaricIsotopeCorrector(method)
            map1 = pyopenms.ConsensusMap()
            map2 = pyopenms.ConsensusMap()
            assert inst.correctIsotopicImpurities(map1, map2)


if __name__ == '__main__':
    unittest.main()
