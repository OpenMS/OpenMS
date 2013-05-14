import unittest
import os

import pdb
import pyopenms

class TestIdXML(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.idXML")

    def test_readfile(self):
        idxml_file = pyopenms.IdXMLFile()
        peps = []
        prots = []
        idxml_file.load(self.filename, prots, peps)


    def test_readfile_content(self):
        idxml_file = pyopenms.IdXMLFile()
        peps = []
        prots = []
        idxml_file.load(self.filename, prots, peps)

        assert len(prots) == 1
        assert len(peps) == 3

if __name__ == '__main__':
    unittest.main()
