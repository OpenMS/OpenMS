from Types cimport *
from String cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set

from MSExperiment cimport *
from FeatureMap cimport *
from Precursor cimport *


cdef extern from "<OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>" namespace "OpenMS":

    cdef cppclass PrecursorCorrection:

        void getPrecursors(MSExperiment & exp,
                           libcpp_vector[ Precursor ] & precursors,
                           libcpp_vector[ double ] & precursors_rt,
                           libcpp_vector[ size_t ] & precursor_scan_index) nogil except +

        void writeHist(String & out_csv,
                       libcpp_vector[ double ] & delta_mzs,
                       libcpp_vector[ double ] & mzs,
                       libcpp_vector[ double ] & rts) nogil except +

        libcpp_set[ size_t ] correctToNearestMS1Peak(MSExperiment & exp,
                                                     double mz_tolerance,
                                                     bool ppm,
                                                     libcpp_vector[ double ] & delta_mzs,
                                                     libcpp_vector[ double ] & mzs,
                                                     libcpp_vector[ double ] & rts) nogil except +

        libcpp_set[ size_t ] correctToHighestIntensityMS1Peak(MSExperiment & exp,
                                                              double mz_tolerance,
                                                              bool ppm,
                                                              libcpp_vector[ double ] & delta_mzs,
                                                              libcpp_vector[ double ] & mzs,
                                                              libcpp_vector[ double ] & rts) nogil except +
        
        libcpp_set[ size_t ] correctToNearestFeature(FeatureMap & features,
                                                     MSExperiment & exp,
                                                     double rt_tolerance_s,
                                                     double mz_tolerance,
                                                     bool ppm,
                                                     bool believe_charge,
                                                     bool keep_original,
                                                     bool all_matching_features,
                                                     int max_trace,
                                                     int debug_level) nogil except +

