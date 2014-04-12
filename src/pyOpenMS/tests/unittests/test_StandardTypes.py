def test_if_available():
    import pyopenms
    assert pyopenms.PeakMap == pyopenms.MSExperiment
    assert pyopenms.RichPeakMap == pyopenms.RichMSExperiment
    assert pyopenms.Chromatorgram == pyopenms.MSChromatogram
    assert pyopenms.PeakSpectrum == pyopenms.MSSpectrum
    assert pyopenms.RichPeakSpectrum == pyopenms.RichMSSpectrum
