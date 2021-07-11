from Types cimport *
from BaseFeature cimport *
from MSExperiment cimport *
from KDTreeFeatureMaps cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector



cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>" namespace "OpenMS":

    cdef cppclass FeatureMapping "OpenMS::FeatureMapping":
       FeatureMapping() nogil except +
       FeatureMapping(FeatureMapping) nogil except +

       FeatureMapping_FeatureToMs2Indices assignMS2IndexToFeature(MSExperiment& spectra,
                                                                  KDTreeFeatureMaps& fp_map_kd,
                                                                  double& precursor_mz_tolerance,
                                                                  double& precursor_rt_tolerance,
                                                                  bool ppm) nogil except +
            # wrap-doc:
            #   Allocate ms2 spectra to feature within the minimal distance
            #   -----
            #   :param spectra: Input of PeakMap/MSExperiment with spectra information
            #   :param fp_map_kd: KDTree used for query and match spectra with features
            #   :param precursor_mz_tolerance: mz_tolerance used for query
            #   :param precursor_rt_tolernace: rt tolerance used for query
            #   :param ppm: mz tolerance window calculation in ppm or Da
            #   :returns: FeatureToMs2Indices

    cdef cppclass FeatureMapping_FeatureToMs2Indices "OpenMS::FeatureMapping::FeatureToMs2Indices":
        FeatureMapping_FeatureToMs2Indices()
        FeatureMapping_FeatureToMs2Indices(FeatureMapping_FeatureToMs2Indices)






