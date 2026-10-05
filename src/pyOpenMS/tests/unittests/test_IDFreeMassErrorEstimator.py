import math

from pyopenms import (
    IDFreeMassErrorEstimator,
    MSSpectrum,
    Precursor,
    SpectrumSettings,
)


def _spectrum(cycle: int, target: int = 0) -> MSSpectrum:
    phase = cycle * 0.73 + target * 0.41
    precursor_mz = 500.0 + target * 20.0 + math.sin(phase) * 0.0007
    fragment_shift = math.sin(phase * 1.31) * 0.0015

    spectrum = MSSpectrum()
    spectrum.setMSLevel(2)
    spectrum.setRT(float(cycle * 10 + target))
    spectrum.setType(SpectrumSettings.SpectrumType.CENTROID)
    spectrum.set_peaks((
        [150.0 + i * 35.0 + fragment_shift for i in range(12)],
        [100.0 + i for i in range(12)],
    ))

    precursor = Precursor()
    precursor.setMZ(precursor_mz)
    precursor.setCharge(2)
    spectrum.setPrecursors([precursor])
    return spectrum


def test_id_free_mass_error_estimator_streaming_binding():
    parameters = IDFreeMassErrorEstimator.Parameters()
    parameters.min_spectrum_pairs = 4
    parameters.min_precursor_clusters = 2
    parameters.min_tolerance_pairs = 4
    parameters.min_tolerance_clusters = 2
    parameters.min_fragment_pairs = 10
    parameters.min_fragment_tolerance_pairs = 10
    parameters.min_fragment_tolerance_spectra = 2

    estimator = IDFreeMassErrorEstimator(parameters)
    for cycle in range(8):
        for target in range(4):
            estimator.consumeSpectrum(_spectrum(cycle, target))

    result = estimator.getResult()
    assert result.precursor_ppm is not None
    assert result.precursor_tolerance_ppm is not None
    assert result.fragment_ppm is not None
    assert result.fragment_tolerance_ppm is not None
    assert result.fragment_tolerance_da is None
    assert (
        result.fragment_resolution_regime
        == IDFreeMassErrorEstimator.FragmentResolutionRegime.HIGH_RESOLUTION
    )
    assert result.precursor_tolerance_ppm.unit == "ppm"
    assert result.fragment_tolerance_ppm.unit == "ppm"
    assert result.diagnostics.precursor_clusters_used == 4
    assert result.diagnostics.fragment_centroid_spectra == 32

    precision = estimator.getPrecursorPrecisionPPM()
    assert precision is not None
    assert precision.single_measurement_sigma > 0.0
