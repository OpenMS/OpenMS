import unittest

import pyopenms


def _feature(rt, mz, intensity, charge):
    f = pyopenms.Feature()
    f.setRT(rt)
    f.setMZ(mz)
    f.setIntensity(intensity)
    f.setCharge(charge)
    return f


def _peptides(sequence):
    """A peptide identification of the search run "search" with one hit"""
    hit = pyopenms.PeptideHit()
    hit.setSequence(pyopenms.AASequence.fromString(sequence))
    pep_id = pyopenms.PeptideIdentification()
    pep_id.setIdentifier("search")
    pep_id.setScoreType("score")
    pep_id.setHits([hit])
    pep_ids = pyopenms.PeptideIdentificationList()
    pep_ids.push_back(pep_id)
    return pep_ids


class TestFeatureOverlapFilter(unittest.TestCase):

    def test_mergeCoincidentFeatures(self):
        fmap = pyopenms.FeatureMap()
        fmap.push_back(_feature(659.206, 698.827042711, 100.0, 2))
        fmap.push_back(_feature(659.206, 698.827042711, 200.0, 2))
        fmap.push_back(_feature(659.206, 698.827042711, 300.0, 3))
        run = pyopenms.ProteinIdentification()
        run.setIdentifier("search")
        ids = pyopenms.MapPeptideIdentifications()
        ids.proteins = [run]
        ids.setFeature([0], _peptides("ESKS(Phospho)SPRPTAEK"))
        ids.setFeature([1], _peptides("ESKSSPRPT(Phospho)AEK"))
        pyopenms.IdentificationDataConverter.importMapIdentifications(fmap, ids)

        self.assertEqual(pyopenms.FeatureOverlapFilter.mergeCoincidentFeatures(fmap), 1)
        self.assertEqual(fmap.size(), 2)
        self.assertAlmostEqual(fmap[0].getIntensity(), 200.0)
        self.assertEqual(len(pyopenms.IdentificationDataConverter.exportMapIdentifications(fmap).at([0])), 2)
        self.assertEqual(fmap[1].getCharge(), 3)


if __name__ == "__main__":
    unittest.main()
