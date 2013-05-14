
import pyopenms
import unittest

class TestLightTargetedExperiment(unittest.TestCase):


    def setUp(self):

        lt = pyopenms.LightTransition()
        lt.charge = 2
        lt.transition_name = "X"
        lt.peptide_ref = "Y"
        lt.library_intensity = 12.0
        lt.product_mz = 22.0
        self.lt = lt

        lm = pyopenms.LightModification()
        lm.location = 13
        lm.unimod_id = "ID"

        self.lm = lm

        lpep = pyopenms.LightPeptide()
        lpep.rt = 12.0
        lpep.charge =  2
        lpep.sequence =  "SEQ"
        lpep.protein_ref =  "REF"


        lpep.modifications = [lm]
        self.lpep = lpep

        lprot = pyopenms.LightProtein()
        lprot.id = "1234"
        lprot.sequence = "ABC"

        self.lprot = lprot

        lte = pyopenms.LightTargetedExperiment()
        lte.peptides = [self.lpep]
        lte.proteins = [self.lprot]
        lte.transitions = [self.lt]

        self.lte = lte


    @staticmethod
    def _test_light_transition(lt):
        assert lt.charge == 2
        assert lt.getProductChargeState() == lt.charge
        assert lt.transition_name == "X"
        assert lt.peptide_ref == "Y"
        assert lt.library_intensity == 12.0
        assert lt.product_mz == 22.0

    @staticmethod
    def _test_light_modification(lm):
        assert lm.location == 13
        assert lm.unimod_id == "ID"

    @staticmethod
    def _test_light_peptide(lpep):

        assert lpep.rt == 12.0
        assert lpep.charge ==  2
        assert lpep.sequence ==  "SEQ"
        assert lpep.protein_ref ==  "REF"

        mod, = lpep.modifications
        TestLightTargetedExperiment._test_light_modification(mod)

    @staticmethod
    def _test_light_protein(lprot):
        assert lprot.id == "1234"
        assert lprot.sequence == "ABC"

    def test_light_transition(self):
        TestLightTargetedExperiment._test_light_transition(self.lt)


    def test_light_modification(self):
        TestLightTargetedExperiment._test_light_modification(self.lm)


    def test_light_peptide(self):
        TestLightTargetedExperiment._test_light_peptide(self.lpep)


    def test_light_protein(self):
        TestLightTargetedExperiment._test_light_protein(self.lprot)

    def test_light_targeted_experiment(self):
        lprot, = self.lte.proteins
        lpep, = self.lte.peptides
        ltrans, = self.lte.transitions

        TestLightTargetedExperiment._test_light_protein(lprot)
        TestLightTargetedExperiment._test_light_peptide(lpep)
        TestLightTargetedExperiment._test_light_transition(ltrans)

        ltrans, = self.lte.getTransitions()
        TestLightTargetedExperiment._test_light_transition(ltrans)

