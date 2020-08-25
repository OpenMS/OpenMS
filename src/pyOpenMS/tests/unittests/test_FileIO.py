import unittest
import os

import pyopenms

class TestPepXML(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.pep.xml").encode()

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

        self.assertEqual( len(prots),  3)
        self.assertEqual( len(peps),  19)

        self.assertEqual( peps[0].getHits()[0].getSequence().toString(), ".(Glu->pyro-Glu)ELNKEMAAEKAKAAAG")

class TestIdXML(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.idXML").encode()

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

        self.assertEqual( len(prots),  1)
        self.assertEqual( len(peps),  3)

class TestIndexedMzMLFileLoader(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.indexed.mzML").encode()

    def test_readfile(self):
        e = pyopenms.OnDiscMSExperiment();
        success = pyopenms.IndexedMzMLFileLoader().load(self.filename, e)

        self.assertTrue(success)

    def test_readfile_content(self):
        e = pyopenms.OnDiscMSExperiment();
        pyopenms.IndexedMzMLFileLoader().load(self.filename, e)

        self.assertEqual( e.getNrSpectra() ,  2)
        self.assertEqual( e.getNrChromatograms() ,  1)

        s = e.getSpectrum(0)
        data_mz, data_int = s.get_peaks()
        self.assertEqual( len(data_mz), 19914)
        self.assertEqual( len(data_int), 19914)

        self.assertEqual( len(e.getSpectrum(1).get_peaks()[0]), 19800)
        self.assertEqual( len(e.getSpectrum(1).get_peaks()[1]), 19800)
        self.assertEqual( len(e.getChromatogram(0).get_peaks()[0]), 48)
        self.assertEqual( len(e.getChromatogram(0).get_peaks()[1]), 48)

        if False:
            # Currently we don't deal with exceptions properly
            raised = False
            try:
                e.getChromatogram(2).get_peaks()
            except Exception as e:
                raised = True

            self.assertTrue(raised)

class TestIndexedMzMLFile(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.indexed.mzML").encode()

    def test_readfile(self):
        f = pyopenms.IndexedMzMLHandler()
        f.openFile(self.filename)

        self.assertTrue(f.getParsingSuccess())

    def test_readfile_content(self):
        f = pyopenms.IndexedMzMLHandler()
        f.openFile(self.filename)

        self.assertEqual( f.getNrSpectra() ,  2)
        self.assertEqual( f.getNrChromatograms() ,  1)

        s = f.getSpectrumById(0)
        mzdata = s.getMZArray()
        intdata = s.getIntensityArray()
        self.assertEqual( len(mzdata), 19914)
        self.assertEqual( len(intdata), 19914)

        s = f.getMSSpectrumById(0)
        self.assertEqual(s.size(), 19914)
        s = f.getMSSpectrumById(1)
        self.assertEqual(s.size(), 19800)
        c = f.getMSChromatogramById(0)
        self.assertEqual(c.size(), 48)

if __name__ == '__main__':
    unittest.main()

