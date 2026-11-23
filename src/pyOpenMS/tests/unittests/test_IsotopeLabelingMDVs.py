"""IsotopeLabelingMDVs updates the Feature passed in.

calculateIsotopicPurity() and calculateMDVAccuracy() used to compute on an empty
Feature created inside the binding. Expected values are those of the C++ class test.
"""
import pytest

import pyopenms as oms


def _feature(peptide_ref, intensities):
    feature = oms.Feature()
    feature.setMetaValue("PeptideRef", peptide_ref)
    subordinates = []
    for i, intensity in enumerate(intensities):
        sub = oms.Feature()
        sub.setMetaValue("native_id", f"{peptide_ref}_{117 + i}")
        sub.setMetaValue("peak_apex_int", intensity)
        subordinates.append(sub)
    feature.setSubordinates(subordinates)
    return feature


def test_calculate_isotopic_purity():
    lactate = _feature("Lactate1", [1.00e+00, 3.324e-05, 2.825e-04, 7.174e-05])
    mdvs = oms.IsotopeLabelingMDVs()
    mdvs.calculateIsotopicPurity(lactate, [0.5, 0.7, 98.8, 0.0, 0.0, 0.0], "1_2-13C_glucose_experiment")
    mdvs.calculateIsotopicPurity(lactate, [0.5, 0.0, 0.1, 0.2, 3.6, 95.5], "U-13C_glucose_experiment")
    assert lactate.getMetaValue("1_2-13C_glucose_experiment") * 100 == pytest.approx(99.6469, rel=1e-4)
    assert lactate.getMetaValue("U-13C_glucose_experiment") * 100 == pytest.approx(99.2517, rel=1e-4)


def test_calculate_mdv_accuracy():
    accoa = _feature("accoa", [0.627, 0.253, 0.096, 0.02, 0.004, 0.001])
    oms.IsotopeLabelingMDVs().calculateMDVAccuracy(accoa, "peak_apex_int", "C23H37N7O17P3S")
    assert accoa.getMetaValue("average_accuracy") == pytest.approx(0.02374, rel=1e-3)
    abs_diff = [sub.getMetaValue("absolute_difference") for sub in accoa.getSubordinates()]
    assert abs_diff == pytest.approx([0.0632108, 0.0505238, 0.0119821, 0.0014131, 0.0000315, 0.0003232],
                                     abs=1e-5)
