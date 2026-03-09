import unittest
import pyopenms as oms


class TestSequenceCoverage(unittest.TestCase):

    def test_partial_coverage(self):
        """Partial sequence coverage return 53.846"""
        protein = oms.AASequence.fromString("PEPTIDEAAAAAA")
        peptides = [
            oms.AASequence.fromString("PEPTIDE")
        ]

        # PEPTIDE covers 7 of 13 amino acids
        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertAlmostEqual(coverage, 53.84615384615385, places=6)

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

    def test_full_coverage_with_repeated_peptide(self):
        """Repeated peptide occurrences cover the full protein"""
        protein = oms.AASequence.fromString("PEPTIDEPEPTIDE")
        peptides = [
            oms.AASequence.fromString("PEPTIDE")
        ]

        coverage = oms.SequenceCoverage.getCoverage(protein, peptides)
        self.assertEqual(coverage, 100.0)


if __name__ == "__main__":
    unittest.main()
