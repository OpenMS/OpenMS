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
        SpectraSTSimilarityScore(SpectraSTSimilarityScore &) nogil except +
        # TODO operator ()
        # double operator()(MSSpectrum & spec1, MSSpectrum & spec2) nogil except +
        # double operator()(BinnedSpectrum & bin1, BinnedSpectrum & bin2) nogil except +
        # double operator()(MSSpectrum & spec) nogil except +
        bool preprocess(MSSpectrum & spec, float remove_peak_intensity_threshold, UInt cut_peaks_below, Size min_peak_number, Size max_peak_number) nogil except +
        # wrap-doc:
                #   Preprocesses the spectrum
                #   -----
                #   The preprocessing removes peak below a intensity threshold, reject spectra that does
                #   not have enough peaks, and cuts peaks exceeding the max_peak_number most intense peaks
                #   -----
                #   :returns: true if spectrum passes filtering

        BinnedSpectrum transform(MSSpectrum & spec) nogil except + # wrap-doc:Spectrum is transformed into a binned spectrum with bin size 1 and spread 1 and the intensities are normalized
        double dot_bias(BinnedSpectrum & bin1, BinnedSpectrum & bin2, double dot_product) nogil except +
        # wrap-doc:
                #   Calculates how much of the dot product is dominated by a few peaks
                #   -----
                #   :param dot_product: If -1 this value will be calculated as well.
                #   :param bin1: First spectrum in binned representation
                #   :param bin2: Second spectrum in binned representation

        double delta_D(double top_hit, double runner_up) nogil except +
        # wrap-doc:
                #   Calculates the normalized distance between top_hit and runner_up
                #   -----
                #   :param top_hit: Is the best score for a given match
                #   :param runner_up: A match with a worse score than top_hit, e.g. the second best score
                #   :returns: normalized distance

        double compute_F(double dot_product, double delta_D, double dot_bias) nogil except +
        # wrap-doc:
                #   Computes the overall all score
                #   -----
                #   :param dot_product: dot_product of a match
                #   :param delta_D: delta_D should be calculated after all dot products for a unidentified spectrum are computed
                #   :param dot_bias
                #   :returns: The SpectraST similarity score

        # POINTER # MSSpectrumCompareFunctor * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Reimplemented from PeakSpectrumCompareFunctor
