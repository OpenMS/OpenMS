from Types cimport *
from libcpp cimport bool
#from PeakSpectrumCompareFunctor cimport *
#from BinnedSpectrum cimport *
from String cimport *
from MSSpectrum cimport *
from Peak1D cimport *

# TODO add BinnedSpectrum
cdef extern from "<OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>" namespace "OpenMS":
    
    cdef cppclass SpectraSTSimilarityScore:
        #  PeakSpectrumCompareFunctor inheritance
        SpectraSTSimilarityScore() nogil except +
        SpectraSTSimilarityScore(SpectraSTSimilarityScore) nogil except +
        # TODO operator ()
        # DoubleReal operator()(MSSpectrum[Peak1D] & spec1, MSSpectrum[Peak1D] & spec2) nogil except +
        # DoubleReal operator()(BinnedSpectrum & bin1, BinnedSpectrum & bin2) nogil except +
        # DoubleReal operator()(MSSpectrum[Peak1D] & spec) nogil except +
        bool preprocess(MSSpectrum[Peak1D] & spec, Real remove_peak_intensity_threshold, UInt cut_peaks_below, Size min_peak_number, Size max_peak_number) nogil except +
        # BinnedSpectrum transform(MSSpectrum[Peak1D] & spec) nogil except +
        # DoubleReal dot_bias(BinnedSpectrum & bin1, BinnedSpectrum & bin2, DoubleReal dot_product) nogil except +
        DoubleReal delta_D(DoubleReal top_hit, DoubleReal runner_up) nogil except +
        DoubleReal compute_F(DoubleReal dot_product, DoubleReal delta_D, DoubleReal dot_bias) nogil except +
        # POINTER # MSSpectrum[Peak1D]CompareFunctor * create() nogil except +
        String getProductName() nogil except +

