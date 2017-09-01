def test_if_available():
    import pyopenms
    assert pyopenms.PeakMap == pyopenms.MSExperiment
    # assert pyopenms.Chromatogram == pyopenms.MSChromatogram
    assert pyopenms.PeakSpectrum == pyopenms.MSSpectrum

