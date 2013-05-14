import unittest
import os

import pyopenms

class TestPepXML(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.pep.xml")

    def test_readfile(self):

        pepxml_file = pyopenms.PepXMLFile()
        peps = []
        prots = []
        pepxml_file.load(self.filename, prots, peps)


    def test_readfile_content(self):

        pepxml_file = pyopenms.PepXMLFile()
        peps = []
        prots = []
        pepxml_file.load(self.filename, prots, peps)

        assert len(prots) == 1
        assert len(peps) == 3

        assert peps[0].getHits()[0].getSequence().toString() == "LAPSAAEDGAFR"

if __name__ == '__main__':
    unittest.main()



