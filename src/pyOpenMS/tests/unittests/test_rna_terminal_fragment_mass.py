"""RNA a1-B ions retain the mass contribution of a 5-prime modification."""

import pyopenms as oms
import pytest


@pytest.mark.parametrize("charge", [-1, 1])
def test_five_prime_a1_base_loss(charge):
    generator = oms.NucleicAcidSpectrumGenerator()
    parameters = generator.getParameters()
    parameters.setValue("add_metainfo", "true")
    parameters.setValue("add_a-B_ions", "true")
    parameters.setValue("add_first_prefix_ion", "true")
    generator.setParameters(parameters)
    sequence = oms.NASequence.fromString("pACGU")
    spectrum = generator.getSpectrum(sequence, charge, charge)
    labels = spectrum.getStringDataArrays()[0]
    indices = [i for i in range(spectrum.size()) if labels[i] == "a1-B"]
    assert len(indices) == 1
    expected = sequence.getPrefix(1).getMonoWeight(oms.NASequence.NASFragmentType.AminusB, charge)
    assert spectrum[indices[0]].getMZ() == pytest.approx(expected)
