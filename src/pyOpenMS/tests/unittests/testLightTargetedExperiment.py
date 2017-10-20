
import pyopenms
import unittest

class TestLightTargetedExperiment(unittest.TestCase):


    def setUp(self):

        lt = pyopenms.LightTransition()
        lt.fragment_charge = 2
        lt.transition_name = b"X"
        lt.peptide_ref = b"Y"
        lt.library_intensity = 12.0
        lt.product_mz = 22.0
        self.lt = lt

        lm = pyopenms.LightModification()
        lm.location = 13
        lm.unimod_id = 4

        self.lm = lm

        lpep = pyopenms.LightCompound()
        lpep.rt = 12.0
        lpep.charge =  2
        lpep.sequence =  b"SEQ"
        lpep.protein_refs =  [b"REF"]


        lpep.modifications = [lm]
        self.lpep = lpep

        lprot = pyopenms.LightProtein()
        lprot.id = b"1234"
        lprot.sequence = b"ABC"

        self.lprot = lprot

        lte = pyopenms.LightTargetedExperiment()
        lte.compounds = [self.lpep]
        lte.proteins = [self.lprot]
        lte.transitions = [self.lt]

        self.lte = lte


    @staticmethod
    def _test_light_transition(lt):
        assert lt.fragment_charge == 2
        assert lt.getProductChargeState() == lt.fragment_charge
        assert lt.transition_name == b"X"
        assert lt.peptide_ref == b"Y"
        assert lt.library_intensity == 12.0
        assert lt.product_mz == 22.0

    @staticmethod
    def _test_light_modification(lm):
        assert lm.location == 13
        assert lm.unimod_id == 4

    @staticmethod
    def _test_light_peptide(lpep):

        assert lpep.rt == 12.0
        assert lpep.charge ==  2
        assert lpep.sequence ==  b"SEQ"
        assert lpep.protein_refs ==  [b"REF"]

        mod, = lpep.modifications
        TestLightTargetedExperiment._test_light_modification(mod)

    @staticmethod
    def _test_light_protein(lprot):
        assert lprot.id == b"1234"
        assert lprot.sequence == b"ABC"

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
        lpep, = self.lte.compounds
        ltrans, = self.lte.transitions

        TestLightTargetedExperiment._test_light_protein(lprot)
        TestLightTargetedExperiment._test_light_peptide(lpep)
        TestLightTargetedExperiment._test_light_transition(ltrans)

        ltrans, = self.lte.getTransitions()
        TestLightTargetedExperiment._test_light_transition(ltrans)

