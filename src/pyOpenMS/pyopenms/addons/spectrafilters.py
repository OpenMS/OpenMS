"""Backward-compatible aliases for spectra filter classes."""
from . import addon


def _make_filter_spectrum():
    """Create filterSpectrum method."""
    def filterSpectrum(self, spectrum):
        """Alias for filterPeakSpectrum (backward compatibility)."""
        return self.filterPeakSpectrum(spectrum)
    return filterSpectrum


# filterPeakSpectrum is also exposed as filterSpectrum for backwards compatibility
for _cls_name in ("NLargest", "Normalizer", "RankScaler", "SqrtScaler",
                   "ThresholdMower", "WindowMower", "GaussFilter",
                   "SavitzkyGolayFilter", "LowessSmoothing",
                   "MorphologicalFilter", "BernNorm"):
    addon(_cls_name)(_make_filter_spectrum())
