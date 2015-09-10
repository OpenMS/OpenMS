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
        # double operator()(MSSpectrum[Peak1D] & spec1, MSSpectrum[Peak1D] & spec2) nogil except +
        # double operator()(BinnedSpectrum & bin1, BinnedSpectrum & bin2) nogil except +
        # double operator()(MSSpectrum[Peak1D] & spec) nogil except +
        bool preprocess(MSSpectrum[Peak1D] & spec, float remove_peak_intensity_threshold, UInt cut_peaks_below, Size min_peak_number, Size max_peak_number) nogil except +
        BinnedSpectrum transform(MSSpectrum[Peak1D] & spec) nogil except +
        double dot_bias(BinnedSpectrum & bin1, BinnedSpectrum & bin2, double dot_product) nogil except +
        double delta_D(double top_hit, double runner_up) nogil except +
        double compute_F(double dot_product, double delta_D, double dot_bias) nogil except +
        # POINTER # MSSpectrum[Peak1D]CompareFunctor * create() nogil except +
        String getProductName() nogil except +

