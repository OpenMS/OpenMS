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

    def test_TargetedExperiment(self):
        targeted = pyopenms.TargetedExperiment();
        tramlfile = pyopenms.TraMLFile();
        tramlfile.load(self.filename, targeted);
        self.assertEqual(len( targeted.getTransitions() ), 3 )

        targeted.setCVs(targeted.getCVs())
        targeted.setTargetCVTerms(targeted.getTargetCVTerms())
        targeted.setPeptides(targeted.getPeptides())
        targeted.setProteins(targeted.getProteins())
        targeted.setTransitions(targeted.getTransitions())

        first_transition = targeted.getTransitions()[0] 
        first_peptide = targeted.getPeptides()[0] 
        targeted.addTransition(first_transition)
        targeted.addPeptide(first_peptide)

        self.assertIsNotNone( targeted.getPeptideByRef(first_transition.getPeptideRef()) )
        self.assertIsNotNone( targeted.getProteinByRef(first_peptide.protein_refs[0]) )


if __name__ == '__main__':
    unittest.main()
