def test_if_available():
    import pyopenms
    assert pyopenms.PeakMap == pyopenms.MSRun
    assert pyopenms.PeakSpectrum == pyopenms.MSSpectrum
