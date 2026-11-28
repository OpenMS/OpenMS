"""
Tests for MSChromatogram class bindings.
"""

import pytest


def test_mschromatogram_copy_meta_data_from():
    """copyMetaDataFrom copies only the meta data; peaks and data arrays stay."""
    from pyopenms import MSChromatogram, ChromatogramSettings, Precursor, FloatDataArray

    source = MSChromatogram()
    source.set_peaks([[1.0], [10.0]])
    source.setName("source")
    source.setNativeID("chrom=7")
    source.setChromatogramType(ChromatogramSettings.ChromatogramType.SELECTED_REACTION_MONITORING_CHROMATOGRAM)
    prec = Precursor()
    prec.setMZ(500.25)
    source.setPrecursor(prec)
    source.setMetaValue("label", "source")
    fda = FloatDataArray()
    fda.push_back(1.0)
    source.setFloatDataArrays([fda])

    target = MSChromatogram()
    target.set_peaks([[5.0, 6.0], [30.0, 40.0]])
    target_fda = FloatDataArray()
    target_fda.setName("target")
    target_fda.push_back(7.0)
    target_fda.push_back(8.0)
    target.setFloatDataArrays([target_fda])

    target.copyMetaDataFrom(source)

    assert target.getName() == "source"
    assert target.getNativeID() == "chrom=7"
    assert target.getChromatogramType() == ChromatogramSettings.ChromatogramType.SELECTED_REACTION_MONITORING_CHROMATOGRAM
    assert target.getPrecursor().getMZ() == pytest.approx(500.25)
    assert target.getMetaValue("label") == "source"

    rt, intensity = target.get_peaks()
    assert list(rt) == pytest.approx([5.0, 6.0])
    assert list(intensity) == pytest.approx([30.0, 40.0])
    target_fdas = target.getFloatDataArrays()
    assert len(target_fdas) == 1
    assert target_fdas[0].getName() == "target"
    assert list(target_fdas[0].get_data()) == pytest.approx([7.0, 8.0])

    assert source.size() == 1
    assert len(source.getFloatDataArrays()) == 1
