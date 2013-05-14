import unittest
import os

import pyopenms

class TestTraMLFile(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.TraML")

    def test_readfile(self):
        targeted = pyopenms.TargetedExperiment();
        tramlfile = pyopenms.TraMLFile();
        tramlfile.load(self.filename, targeted);

    def test_readfile_content(self):
        targeted = pyopenms.TargetedExperiment();
        tramlfile = pyopenms.TraMLFile();
        tramlfile.load(self.filename, targeted);
        self.assertEqual(len( targeted.getTransitions() ), 3 )

        self.assertAlmostEqual(targeted.getTransitions()[0].getPrecursorMZ(), 500.0)
        self.assertAlmostEqual(targeted.getTransitions()[0].getProductMZ(), 628.45, places=4)
        self.assertEqual(targeted.getTransitions()[0].getName(), "tr1" )
        self.assertEqual(targeted.getTransitions()[0].getNativeID(), "tr1" )
        self.assertEqual(targeted.getTransitions()[0].getPeptideRef(), "tr_gr1")

if __name__ == '__main__':
    unittest.main()
