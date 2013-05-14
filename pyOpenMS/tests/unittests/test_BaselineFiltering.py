import unittest
import os

import pyopenms

class TestMorphologicalFilter(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML")
        self.exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, self.exp)

    def test_init(self):
        thisfilter = pyopenms.MorphologicalFilter();

    def test_run(self):
        thisfilter = pyopenms.MorphologicalFilter();
        old_firstspec = self.exp[0]
        # needs different parameters to have any effect ...
        params = pyopenms.MorphologicalFilter().getDefaults();
        params.setValue("struc_elem_length", pyopenms.DataValue(0.05), '')
        thisfilter.setParameters(params);
        thisfilter.filterExperiment(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # MZ should not change, Intensity should
        self.assertEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

if __name__ == '__main__':
    unittest.main()
