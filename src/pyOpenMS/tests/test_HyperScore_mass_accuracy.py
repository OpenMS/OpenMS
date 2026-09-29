"""Mass-accuracy HyperScore binding: score and unweighted detail semantics."""

import numpy as np
import pytest
import pyopenms as oms


def test_mass_accuracy_score_and_detail():
    generator = oms.TheoreticalSpectrumGenerator()
    params = generator.getParameters()
    params.setValue("add_metainfo", "true")
    generator.setParameters(params)
    theoretical = oms.MSSpectrum()
    generator.getSpectrum(theoretical, oms.AASequence.fromString("PEPTIDEK"), 1, 1)
    exact = oms.HyperScore.computeMassAccuracy(20.0, True, theoretical, theoretical, 7.0)
    assert len(exact) == 4
    assert exact[0] == pytest.approx(oms.HyperScore.compute(20.0, True, theoretical, theoretical))
    assert exact[1] > 0 and exact[2] > 0
    assert exact[3] == pytest.approx(0.0)

    observed = oms.MSSpectrum()
    mz, intensity = theoretical.get_peaks()
    observed.set_peaks((mz * (1.0 + 7e-6), intensity))
    shifted = oms.HyperScore.computeMassAccuracy(20.0, True, observed, theoretical, 7.0)
    assert 0 < shifted[0] < exact[0]
    assert shifted[1:3] == exact[1:3]
    assert shifted[3] == pytest.approx(7.0)
    for invalid in (0.0, -1.0, np.inf, np.nan):
        with pytest.raises(Exception):
            oms.HyperScore.computeMassAccuracy(20.0, True, observed, theoretical, invalid)
