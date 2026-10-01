"""FragmentIonLikelihoodModel binding: training, scoring and the static helpers."""

import pytest
import pyopenms as oms


def _theoretical(sequence):
    generator = oms.TheoreticalSpectrumGenerator()
    params = generator.getParameters()
    params.setValue("add_metainfo", "true")
    generator.setParameters(params)
    spectrum = oms.MSSpectrum()
    generator.getSpectrum(spectrum, oms.AASequence.fromString(sequence), 1, 1)
    return spectrum


def test_train_and_score():
    peptide = "DFPIANGER"
    theoretical = _theoretical(peptide)
    reversed_theoretical = _theoretical("EGNAIPFDR")
    ranks = oms.FragmentIonLikelihoodModel.intensityRanks(theoretical)
    assert len(ranks) == theoretical.size()
    assert min(ranks) == 1 and max(ranks) == theoretical.size()

    model = oms.FragmentIonLikelihoodModel()
    assert not model.isTrained()
    assert model.pseudoCount() == pytest.approx(20.0)
    model.addObservations(theoretical, ranks, theoretical, len(peptide), 2, 20.0, True, False)
    model.addObservations(theoretical, ranks, reversed_theoretical, len(peptide), 2, 20.0, True, True)
    assert model.signalPsms() == 1 and model.noisePsms() == 1
    with pytest.raises(Exception):
        model.score(theoretical, ranks, theoretical, len(peptide), 2, 20.0, True)
    model.finalize()
    assert model.isTrained()

    features = model.score(theoretical, ranks, theoretical, len(peptide), 2, 20.0, True)
    assert features.theoretical_ions == features.matched_ions > 0
    assert features.log_likelihood_ratio > 0.0
    assert features.explained_presence == pytest.approx(1.0)
    assert features.top_predicted_observed == pytest.approx(1.0)
    noise = model.score(theoretical, ranks, reversed_theoretical, len(peptide), 2, 20.0, True, 4)
    assert noise.matched_ions < features.matched_ions
    assert noise.log_likelihood_ratio < features.log_likelihood_ratio

    context = oms.FragmentIonLikelihoodModel.contextOf(False, 2, 1, 3, len(peptide))
    assert context.series == 1 and context.position_bin == 3
    assert 0.0 < model.presenceProbability(context) <= 1.0
    absent = oms.FragmentIonLikelihoodModel.ABSENT
    assert model.logLikelihoodRatio(context, 0) > model.logLikelihoodRatio(context, absent)

    # Contexts outside the model's tables are rejected instead of read.
    invalid = oms.FragmentIonLikelihoodModel.contextOf(False, 2, 1, 3, len(peptide))
    invalid.series = 2
    with pytest.raises(Exception):
        model.presenceProbability(invalid)
    with pytest.raises(Exception):
        model.logLikelihoodRatio(invalid, 0)


def test_static_helpers():
    assert oms.FragmentIonLikelihoodModel.parseIonName("y3++") == (True, False, 3)
    assert oms.FragmentIonLikelihoodModel.parseIonName("b5+") == (True, True, 5)
    assert oms.FragmentIonLikelihoodModel.parseIonName("[M+H]+")[0] is False
    assert oms.FragmentIonLikelihoodModel.rankOutcome(1) == 0
    assert oms.FragmentIonLikelihoodModel.rankOutcome(1000) < oms.FragmentIonLikelihoodModel.ABSENT
    assert oms.FragmentIonLikelihoodModel.OUTCOMES == oms.FragmentIonLikelihoodModel.ABSENT + 1
    assert oms.FragmentIonLikelihoodModel(5.0).pseudoCount() == pytest.approx(5.0)
    with pytest.raises(Exception):
        oms.FragmentIonLikelihoodModel(0.0)
