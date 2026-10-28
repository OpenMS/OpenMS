import unittest

import pyopenms


def _feature(rt, mz, intensity, charge, sequence=None):
    f = pyopenms.Feature()
    f.setRT(rt)
    f.setMZ(mz)
    f.setIntensity(intensity)
    f.setCharge(charge)
    if sequence is not None:
        hit = pyopenms.PeptideHit()
        hit.setSequence(pyopenms.AASequence.fromString(sequence))
        pep_id = pyopenms.PeptideIdentification()
        pep_id.setHits([hit])
        pep_ids = pyopenms.PeptideIdentificationList()
        pep_ids.push_back(pep_id)
        f.setPeptideIdentifications(pep_ids)
    return f


class TestFeatureOverlapFilter(unittest.TestCase):

    def test_mergeCoincidentFeatures(self):
        fmap = pyopenms.FeatureMap()
        fmap.push_back(_feature(659.206, 698.827042711, 100.0, 2, "ESKS(Phospho)SPRPTAEK"))
        fmap.push_back(_feature(659.206, 698.827042711, 200.0, 2, "ESKSSPRPT(Phospho)AEK"))
        fmap.push_back(_feature(659.206, 698.827042711, 300.0, 3))

        self.assertEqual(pyopenms.FeatureOverlapFilter.mergeCoincidentFeatures(fmap), 1)
        self.assertEqual(fmap.size(), 2)
        self.assertAlmostEqual(fmap[0].getIntensity(), 200.0)
        self.assertEqual(len(fmap[0].getPeptideIdentifications()), 2)
        self.assertEqual(fmap[1].getCharge(), 3)


if __name__ == "__main__":
    unittest.main()
