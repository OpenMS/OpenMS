import unittest
import pyopenms as oms


class TestSequenceCoverage(unittest.TestCase):

    def test_peptide_multiple_occurrences(self):
        """Test peptide mapping multiple times covers all occurrences"""
        protein = oms.AASequence.fromString("PEPTIDEPEPTIDE")
        peptides = [
            oms.AASequence.fromString("PEPTIDE")
        ]
        
        # PEPTIDE occurs twice and covers the entire protein
        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertAlmostEqual(coverage, 100.0, places=6)

    def test_partial_coverage(self):
        """Partial sequence coverage is computed correctly"""
        protein = oms.AASequence.fromString("PEPTIDEAAAAAAA")
        peptides = [
            oms.AASequence.fromString("PEPTIDE")
        ]

        # PEPTIDE covers 7 of 14 amino acids (single occurrence)
        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertAlmostEqual(coverage, 50.0, places=6)

    def test_empty_peptide_list(self):
        """Empty peptide list returns zero coverage"""
        protein = oms.AASequence.fromString("PEPTIDE")
        peptides = []

        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertEqual(coverage, 0.0)

    def test_empty_protein(self):
        """Empty protein returns zero coverage"""
        protein = oms.AASequence.fromString("")
        peptides = [oms.AASequence.fromString("PEPTIDE")]

        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertEqual(coverage, 0.0)

    def test_full_coverage(self):
        """Full peptide coverage returns 100 percent"""
        protein = oms.AASequence.fromString("PEPTIDE")
        peptides = [oms.AASequence.fromString("PEPTIDE")]

        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertEqual(coverage, 100.0)


if __name__ == "__main__":
    unittest.main()
