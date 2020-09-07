from Types cimport *
from libcpp cimport bool
#from PeakSpectrumCompareFunctor cimport *
from BinnedSpectrum cimport *
from String cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>" namespace "OpenMS":
    
    cdef cppclass SpectraSTSimilarityScore:
        #  PeakSpectrumCompareFunctor inheritance
        SpectraSTSimilarityScore() nogil except +
        SpectraSTSimilarityScore(SpectraSTSimilarityScore) nogil except +
        # TODO operator ()
        # double operator()(MSSpectrum & spec1, MSSpectrum & spec2) nogil except +
        # double operator()(BinnedSpectrum & bin1, BinnedSpectrum & bin2) nogil except +
        # double operator()(MSSpectrum & spec) nogil except +
        bool preprocess(MSSpectrum & spec, float remove_peak_intensity_threshold, UInt cut_peaks_below, Size min_peak_number, Size max_peak_number) nogil except +
        BinnedSpectrum transform(MSSpectrum & spec) nogil except +
        double dot_bias(BinnedSpectrum & bin1, BinnedSpectrum & bin2, double dot_product) nogil except +
        double delta_D(double top_hit, double runner_up) nogil except +
        double compute_F(double dot_product, double delta_D, double dot_bias) nogil except +
        # POINTER # MSSpectrumCompareFunctor * create() nogil except +
        String getProductName() nogil except +

