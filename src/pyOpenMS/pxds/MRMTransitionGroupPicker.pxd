from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MRMTransitionGroup cimport *
from MRMFeature cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from ChromatogramPeak cimport *
from SignalToNoiseEstimatorMedian cimport *
from SavitzkyGolayFilter cimport *
from GaussFilter cimport *
from LinearResampler cimport *
from PeakPickerHiRes cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>" namespace "OpenMS":
    
    cdef cppclass MRMTransitionGroupPicker(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MRMTransitionGroupPicker() except + nogil 
        MRMTransitionGroupPicker(MRMTransitionGroupPicker &) except + nogil 

        void pickTransitionGroup(MRMTransitionGroup[MSChromatogram, LightTransition] transition_group) except + nogil 
        void pickTransitionGroup(MRMTransitionGroup[MSChromatogram, ReactionMonitoringTransition] transition_group) except + nogil 

        MRMFeature createMRMFeature(MRMTransitionGroup[ MSChromatogram, LightTransition] transition_group,
                                    libcpp_vector[ MSChromatogram ] & picked_chroms,
                                    libcpp_vector[ MSChromatogram ] & smoothed_chroms,
                                    const int chr_idx,
                                    const int peak_idx) except + nogil 

        void remove_overlapping_features(libcpp_vector[ MSChromatogram ] & picked_chroms,
                                         double best_left, 
                                         double best_right) except + nogil 

        void findLargestPeak(libcpp_vector[ MSChromatogram ] & picked_chroms, int & chr_idx, int & peak_idx) except + nogil 

        void findWidestPeakIndices(libcpp_vector[ MSChromatogram ] & picked_chroms, Int & chrom_idx, Int & point_idx) except + nogil 
