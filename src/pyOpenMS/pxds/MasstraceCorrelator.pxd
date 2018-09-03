from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from ConsensusMap cimport *
from MSExperiment cimport *

# typedef std::vector<std::pair<double, double> > MasstracePointsType;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>" namespace "OpenMS":
    
    ctypedef libcpp_vector[libcpp_pair[double,double ] ] MasstracePointsType

    cdef cppclass MasstraceCorrelator(DefaultParamHandler,ProgressLogger) :
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        MasstraceCorrelator() nogil except +
        MasstraceCorrelator(MasstraceCorrelator) nogil except + #wrap-ignore
        void createPseudoSpectra(const ConsensusMap & map_,
                                 MSExperiment & pseudo_spectra,
                                 Size min_peak_nr,
                                 double min_correlation,
                                 int max_lag,
                                 double max_rt_apex_difference) nogil except +

        # void scoreHullpoints(const MasstracePointsType & hull_points1,
        #                      const MasstracePointsType & hull_points2,
        #                      int lag,
        #                      double lag_intensity,
        #                      double pearson_score,
        #                      double min_corr,
        #                      int max_lag,
        #                      double mindiff) nogil except +

        # void createConsensusMapCache(const ConsensusMap & map_,
        #                              libcpp_vector[ MasstracePointsType ] & feature_points,
        #                              libcpp_vector[ libcpp_pair[ double, double ] ] & max_intensities,
        #                              libcpp_vector[ double ] & rt_cache) nogil except +

