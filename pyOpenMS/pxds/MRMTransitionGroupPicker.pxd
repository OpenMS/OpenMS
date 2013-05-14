from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MRMTransitionGroup cimport *
from MRMFeature cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from ChromatogramPeak cimport *
from SignalToNoiseEstimatorMedian cimport *
from SavitzkyGolayFilter cimport *
from GaussFilter cimport *
from LinearResampler cimport *
# from LinearResamplerAlign cimport *
from PeakPickerHiRes cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>" namespace "OpenMS":
    
    cdef cppclass MRMTransitionGroupPicker(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MRMTransitionGroupPicker() nogil except +
        MRMTransitionGroupPicker(MRMTransitionGroupPicker) nogil except + #wrap-ignore
        # TEMPLATE # void pickTransitionGroup(MRMTransitionGroup[ SpectrumT, TransitionT ] & transition_group) nogil except +
        void pickChromatogram(MSSpectrum[ChromatogramPeak] & chromatogram, MSSpectrum[ChromatogramPeak] & smoothed_chrom, MSSpectrum[ChromatogramPeak] & picked_chrom) nogil except +
        # TEMPLATE # MRMFeature createMRMFeature(MRMTransitionGroup[ SpectrumT, TransitionT ] & transition_group, libcpp_vector[ SpectrumT ] & picked_chroms, int & chr_idx, int & peak_idx) nogil except +
        # TEMPLATE # void remove_overlapping_features(libcpp_vector[ SpectrumT ] & picked_chroms, double best_left, double best_right) nogil except +
        void findLargestPeak(libcpp_vector[ MSSpectrum[ChromatogramPeak] ] & picked_chroms, int & chr_idx, int & peak_idx) nogil except +

