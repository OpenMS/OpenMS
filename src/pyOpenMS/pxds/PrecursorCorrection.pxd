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

        String csv_header;

        void getPrecursor(MSExperiment & exp,
                          libcpp_vector[ Precursor ] & precursors,
                          libcpp_vector[ double ] & precursors_rt,
                          libcpp_vector[ Size ] & precursor_scan_index) nogil+

        void writeHist(String & out_csv,
                       libcpp_vector[ double ] & delta_mzs,
                       libcpp_vector[ double ] & mzs,
                       libcpp_vector[ double ] & rts) nogil+

        libcpp_set[ Size ] correctToNearestMS1Peak(MSExperiment & exp,
                                                   double mz_tolerance,
                                                   bool ppm,
                                                   libcpp_vector[ double ] & delta_mzs,
                                                   libcpp_vector[ double ] & mzs,
                                                   libcpp_vector[ double ] & rts) nogil+

        libcpp_set[ Size ] correctToHighestintensityMS1Peak(MSExperiment & exp,
                                                            double mz_tolerance,
                                                            libcpp_vector[ double ] & delta_mzs,
                                                            libcpp_vector[ double ] & mzs,
                                                            libcpp_vector[ double ] & rts) nogil+

        libcpp_set[ Size ] corrrectToNearestFeature(FeatureMap & features,
                                                    MSExperiment & exp,
                                                    double rt_tolerance_s = 0.0,
                                                    double mz_tolerance = 0.0,
                                                    bool ppm = true,
                                                    bool believe_charge = false,
                                                    bool keep_original = false,
                                                    bool all_matching_features = false,
                                                    int max_trace = 2,
                                                    int debug_level = 0) nogil+

