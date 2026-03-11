import unittest
import pyopenms as oms

class TestIDMapper(unittest.TestCase):

    def create_test_feature(self, rt, mz, intensity=1000.0, charge=2):
        f = oms.Feature()
        f.setRT(rt)
        f.setMZ(mz)
        f.setIntensity(intensity)
        f.setCharge(charge)
        return f

    def create_test_peptide_id(self, rt, mz, sequence="PEPTIDER", score=100.0, charge=2):
        pi = oms.PeptideIdentification()
        pi.setRT(rt)
        pi.setMZ(mz)
        pi.setScoreType("search_engine_score")
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(sequence))
        hit.setScore(score)
        hit.setCharge(charge)
        pi.setHits([hit])
        return pi

    def create_test_spectrum(self, rt, ms_level, precursor_mz, precursor_charge=2):
        s = oms.MSSpectrum()
        s.setRT(rt)
        s.setMSLevel(ms_level)
        p = oms.Precursor()
        p.setMZ(precursor_mz)
        p.setCharge(precursor_charge)
        s.setPrecursors([p])
        return s

    def assert_feature_has_identification(self, feature, expected_sequence):
        pep_ids = feature.getPeptideIdentifications()
        self.assertGreater(len(pep_ids), 0, "Feature has no peptide identifications.")
        hits = pep_ids[0].getHits()
        self.assertGreater(len(hits), 0, "Peptide identification has no hits.")
        self.assertEqual(hits[0].getSequence().toString(), expected_sequence)

    def to_pep_list(self, peps):
        """Helper to convert Python lists to C++ PeptideIdentificationList."""
        pep_list = oms.PeptideIdentificationList()
        for p in peps:
            pep_list.append(p)
        return pep_list

    def test_annotate_featuremap_empty_msexperiment(self):
        """Test annotate with an empty MSExperiment."""
        mapper = oms.IDMapper()

        params = mapper.getParameters()
        params.setValue("rt_tolerance", 5.0)
        params.setValue("mz_tolerance", 0.5)
        params.setValue("mz_measure", "Da")
        mapper.setParameters(params)

        features = oms.FeatureMap()
        features.push_back(self.create_test_feature(500.0, 800.0))

        peptide_ids = self.to_pep_list([self.create_test_peptide_id(500.0, 800.0, "TESTPEPTIDE")])

        empty_exp = oms.MSExperiment()
        empty_exp.updateRanges()

        mapper.annotate(features, peptide_ids, [], True, True, empty_exp)
        self.assert_feature_has_identification(features[0], "TESTPEPTIDE")

    def test_annotate_featuremap_with_msexperiment(self):
        """Test with populated MSExperiment containing MS2 spectra."""
        mapper = oms.IDMapper()

        params = mapper.getParameters()
        params.setValue("rt_tolerance", 5.0)
        params.setValue("mz_tolerance", 0.5)
        params.setValue("mz_measure", "Da")
        mapper.setParameters(params)

        features = oms.FeatureMap()
        features.push_back(self.create_test_feature(500.0, 800.0))

        peptide_ids = self.to_pep_list([self.create_test_peptide_id(500.0, 800.0, "MATCHEDPEP")])

        exp = oms.MSExperiment()
        exp.addSpectrum(self.create_test_spectrum(500.0, 2, 800.0))
        exp.updateRanges()

        mapper.annotate(features, peptide_ids, [], True, True, exp)
        self.assert_feature_has_identification(features[0], "MATCHEDPEP")

    def test_annotate_tolerance_settings(self):
        """Test RT and m/z tolerance boundaries."""
        mapper = oms.IDMapper()

        params = mapper.getParameters()
        params.setValue("rt_tolerance", 5.0)
        params.setValue("mz_tolerance", 0.5)
        params.setValue("mz_measure", "Da")
        mapper.setParameters(params)

        features = oms.FeatureMap()
        features.push_back(self.create_test_feature(500.0, 800.0))

        peptide_ids = self.to_pep_list([
            self.create_test_peptide_id(504.0, 800.2, "INSIDE"),
            self.create_test_peptide_id(506.0, 800.6, "OUTSIDE")
        ])

        empty_exp = oms.MSExperiment()
        empty_exp.updateRanges()
        mapper.annotate(features, peptide_ids, [], True, True, empty_exp)

        pep_ids = features[0].getPeptideIdentifications()
        self.assertEqual(len(pep_ids), 1)
        self.assertEqual(pep_ids[0].getHits()[0].getSequence().toString(), "INSIDE")

    def test_annotate_centroid_parameters(self):
        """Test boolean parameter combinations for use_centroid_rt and use_centroid_mz."""
        mapper = oms.IDMapper()

        params = mapper.getParameters()
        params.setValue("rt_tolerance", 5.0)
        params.setValue("mz_tolerance", 0.5)
        params.setValue("mz_measure", "Da")
        mapper.setParameters(params)

        f = self.create_test_feature(500.0, 800.0)

        hull = oms.ConvexHull2D()
        hull.addPoint([490.0, 799.0])
        hull.addPoint([510.0, 799.0])
        hull.addPoint([510.0, 801.0])
        hull.addPoint([490.0, 801.0])

        hulls = f.getConvexHulls()
        hulls.append(hull)
        f.setConvexHulls(hulls)

        combinations = [(True, True), (True, False), (False, True), (False, False)]

        empty_exp = oms.MSExperiment()
        empty_exp.updateRanges()

        for use_centroid_rt, use_centroid_mz in combinations:
            with self.subTest(use_centroid_rt=use_centroid_rt, use_centroid_mz=use_centroid_mz):
                features = oms.FeatureMap()
                features.push_back(f)

                peptide_ids = self.to_pep_list([self.create_test_peptide_id(500.0, 800.0, "HULLTEST")])

                mapper.annotate(features, peptide_ids, [], use_centroid_rt, use_centroid_mz, empty_exp)

                pep_ids = features[0].getPeptideIdentifications()
                self.assertEqual(len(pep_ids), 1)
                self.assertEqual(pep_ids[0].getHits()[0].getSequence().toString(), "HULLTEST")

    def test_annotate_empty_featuremap(self):
        """Verify no exceptions on empty FeatureMap."""
        mapper = oms.IDMapper()
        features = oms.FeatureMap()
        peptide_ids = self.to_pep_list([self.create_test_peptide_id(500.0, 800.0)])

        empty_exp = oms.MSExperiment()
        empty_exp.updateRanges()

        mapper.annotate(features, peptide_ids, [], True, True, empty_exp)

    def test_annotate_empty_peptide_ids(self):
        """Verify no exceptions on empty PeptideIdentifications."""
        mapper = oms.IDMapper()
        features = oms.FeatureMap()
        features.push_back(self.create_test_feature(500.0, 800.0))
        peptide_ids = oms.PeptideIdentificationList()

        empty_exp = oms.MSExperiment()
        empty_exp.updateRanges()

        mapper.annotate(features, peptide_ids, [], True, True, empty_exp)
        self.assertEqual(len(features[0].getPeptideIdentifications()), 0)

    def test_annotate_multiple_features_multiple_ids(self):
        """Test many-to-many mapping behavior."""
        mapper = oms.IDMapper()

        params = mapper.getParameters()
        params.setValue("rt_tolerance", 5.0)
        params.setValue("mz_tolerance", 0.5)
        params.setValue("mz_measure", "Da")
        mapper.setParameters(params)

        features = oms.FeatureMap()
        features.push_back(self.create_test_feature(100.0, 400.0))
        features.push_back(self.create_test_feature(200.0, 500.0))

        peptide_ids = self.to_pep_list([
            self.create_test_peptide_id(100.0, 400.0, "PEPONE"),
            self.create_test_peptide_id(200.0, 500.0, "PEPTWO")
        ])

        empty_exp = oms.MSExperiment()
        empty_exp.updateRanges()

        mapper.annotate(features, peptide_ids, [], True, True, empty_exp)

        self.assert_feature_has_identification(features[0], "PEPONE")
        self.assert_feature_has_identification(features[1], "PEPTWO")

if __name__ == '__main__':
    unittest.main()