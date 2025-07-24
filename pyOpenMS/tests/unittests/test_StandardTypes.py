def test_if_available():
    import pyopenms
    assert pyopenms.PeakMap == pyopenms.MSExperiment
    assert pyopenms.PeakSpectrum == pyopenms.MSSpectrum
