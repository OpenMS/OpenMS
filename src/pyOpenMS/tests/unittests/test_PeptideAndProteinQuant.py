"""
Tests for PeptideAndProteinQuant class methods getPeptideResults() and getProteinResults().

These tests verify that the peptide and protein quantification results are correctly
exposed to Python, addressing GitHub issue #5657.
"""
import unittest
import os

import pyopenms


class TestPeptideAndProteinQuant(unittest.TestCase):
    """Test PeptideAndProteinQuant methods."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.dirname = os.path.dirname(os.path.abspath(__file__))
        cls.feature_file = os.path.join(cls.dirname, "ProteinQuantifier_input.featureXML")
        cls.consensus_file = os.path.join(cls.dirname, "ProteinQuantifier_input.consensusXML")

    def test_init(self):
        """Test PeptideAndProteinQuant can be instantiated."""
        quant = pyopenms.PeptideAndProteinQuant()
        self.assertIsNotNone(quant)

    def test_get_peptide_results_from_feature_map(self):
        """Test getPeptideResults() returns dict with correct peptide data from FeatureMap."""
        # Load feature map
        features = pyopenms.FeatureMap()
        pyopenms.FeatureXMLFile().load(self.feature_file, features)

        # Create experimental design from feature map
        design = pyopenms.ExperimentalDesign.fromFeatureMap(features)

        # Set up quantifier with include_all to get all results
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(features, design)
        quant.quantifyPeptides()

        # Get peptide results
        peptide_results = quant.getPeptideResults()

        # Verify results
        self.assertIsInstance(peptide_results, dict)
        self.assertGreater(len(peptide_results), 0)

        # From C++ test: pep_quant.size() == 7
        self.assertEqual(len(peptide_results), 7)

        # Convert keys to string for easier lookup
        peptide_by_seq = {k.toString(): v for k, v in peptide_results.items()}

        # Check specific peptide AAAAA
        self.assertIn("AAAAA", peptide_by_seq)
        aaaaa_data = peptide_by_seq["AAAAA"]
        self.assertIsNotNone(aaaaa_data)
        # From C++ test: pep_data.psm_count == 2
        self.assertEqual(aaaaa_data.psm_count, 2)
        # From C++ test: pep_data.accessions.size() == 1
        accessions = aaaaa_data.accessions
        self.assertEqual(len(accessions), 1)

        # Check specific peptide CCCCC
        self.assertIn("CCCCC", peptide_by_seq)
        ccccc_data = peptide_by_seq["CCCCC"]
        self.assertEqual(ccccc_data.psm_count, 2)

    def test_get_protein_results_from_feature_map(self):
        """Test getProteinResults() returns dict with correct protein data from FeatureMap."""
        # Load feature map
        features = pyopenms.FeatureMap()
        pyopenms.FeatureXMLFile().load(self.feature_file, features)

        # Create experimental design from feature map
        design = pyopenms.ExperimentalDesign.fromFeatureMap(features)

        # Set up quantifier
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(features, design)
        quant.quantifyPeptides()
        quant.quantifyProteins()

        # Get protein results
        protein_results = quant.getProteinResults()

        # Verify results
        self.assertIsInstance(protein_results, dict)
        self.assertGreater(len(protein_results), 0)

        # From C++ test: prot_quant.size() == 2
        self.assertEqual(len(protein_results), 2)

        # Check specific protein Protein0
        self.assertIn("Protein0", protein_results)
        protein0_data = protein_results["Protein0"]
        self.assertIsNotNone(protein0_data)
        # From C++ test: prot_data.psm_count == 6
        self.assertEqual(protein0_data.psm_count, 6)

        # Check specific protein Protein1
        self.assertIn("Protein1", protein_results)
        protein1_data = protein_results["Protein1"]
        # From C++ test: prot_data.psm_count == 2
        self.assertEqual(protein1_data.psm_count, 2)

    def test_get_peptide_results_from_consensus_map(self):
        """Test getPeptideResults() returns dict with correct peptide data from ConsensusMap."""
        # Load consensus map
        consensus = pyopenms.ConsensusMap()
        pyopenms.ConsensusXMLFile().load(self.consensus_file, consensus)

        # Create experimental design from consensus map
        design = pyopenms.ExperimentalDesign.fromConsensusMap(consensus)

        # Set up quantifier
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(consensus, design)
        quant.quantifyPeptides()

        # Get peptide results
        peptide_results = quant.getPeptideResults()

        # Verify results
        self.assertIsInstance(peptide_results, dict)
        # From C++ test: pep_quant.size() == 4
        self.assertEqual(len(peptide_results), 4)

        # Convert keys to string for easier lookup
        peptide_by_seq = {k.toString(): v for k, v in peptide_results.items()}

        # Check specific peptide AAAK
        self.assertIn("AAAK", peptide_by_seq)
        aaak_data = peptide_by_seq["AAAK"]
        self.assertEqual(aaak_data.psm_count, 1)

    def test_get_protein_results_from_consensus_map(self):
        """Test getProteinResults() returns dict with correct protein data from ConsensusMap."""
        # Load consensus map
        consensus = pyopenms.ConsensusMap()
        pyopenms.ConsensusXMLFile().load(self.consensus_file, consensus)

        # Create experimental design from consensus map
        design = pyopenms.ExperimentalDesign.fromConsensusMap(consensus)

        # Set up quantifier
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(consensus, design)
        quant.quantifyPeptides()
        quant.quantifyProteins()

        # Get protein results
        protein_results = quant.getProteinResults()

        # Verify results
        self.assertIsInstance(protein_results, dict)
        # From C++ test: prot_quant.size() == 1
        self.assertEqual(len(protein_results), 1)

        # Check specific protein
        self.assertIn("Protein", protein_results)
        protein_data = protein_results["Protein"]
        # From C++ test: prot_data.psm_count == 4
        self.assertEqual(protein_data.psm_count, 4)

    def test_get_statistics(self):
        """Test getStatistics() returns correct statistics."""
        # Load feature map
        features = pyopenms.FeatureMap()
        pyopenms.FeatureXMLFile().load(self.feature_file, features)

        # Create experimental design from feature map
        design = pyopenms.ExperimentalDesign.fromFeatureMap(features)

        # Set up quantifier
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(features, design)
        quant.quantifyPeptides()
        quant.quantifyProteins()

        # Get statistics
        stats = quant.getStatistics()

        # Verify statistics (from C++ test)
        self.assertEqual(stats.n_samples, 1)
        self.assertEqual(stats.quant_proteins, 2)
        self.assertEqual(stats.too_few_peptides, 1)
        self.assertEqual(stats.quant_peptides, 5)
        self.assertEqual(stats.total_peptides, 7)
        self.assertEqual(stats.quant_features, 7)
        self.assertEqual(stats.total_features, 8)
        self.assertEqual(stats.blank_features, 0)
        self.assertEqual(stats.ambig_features, 1)

    def test_peptide_data_accessions(self):
        """Test that peptide data accessions are properly exposed as a set."""
        # Load feature map
        features = pyopenms.FeatureMap()
        pyopenms.FeatureXMLFile().load(self.feature_file, features)

        # Create experimental design from feature map
        design = pyopenms.ExperimentalDesign.fromFeatureMap(features)

        # Set up quantifier
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(features, design)
        quant.quantifyPeptides()

        # Get peptide results
        peptide_results = quant.getPeptideResults()

        # Convert keys to string for easier lookup
        peptide_by_seq = {k.toString(): v for k, v in peptide_results.items()}

        # Check GGGGG which maps to 2 proteins (from C++ test)
        self.assertIn("GGGGG", peptide_by_seq)
        ggggg_data = peptide_by_seq["GGGGG"]
        accessions = ggggg_data.accessions
        self.assertEqual(len(accessions), 2)

    def test_peptide_results_key_is_aasequence(self):
        """Test that getPeptideResults() returns AASequence objects as keys."""
        # Load feature map
        features = pyopenms.FeatureMap()
        pyopenms.FeatureXMLFile().load(self.feature_file, features)

        # Create experimental design from feature map
        design = pyopenms.ExperimentalDesign.fromFeatureMap(features)

        # Set up quantifier
        quant = pyopenms.PeptideAndProteinQuant()
        params = quant.getParameters()
        params.setValue("top:include_all", "true")
        quant.setParameters(params)

        # Process data
        quant.readQuantData(features, design)
        quant.quantifyPeptides()

        # Get peptide results
        peptide_results = quant.getPeptideResults()

        # Verify keys are AASequence objects
        for key in peptide_results.keys():
            self.assertIsInstance(key, pyopenms.AASequence)
            # Verify AASequence methods work
            self.assertIsNotNone(key.toString())
            self.assertGreater(key.size(), 0)


if __name__ == '__main__':
    unittest.main()
