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

    cdef cppclass FeatureMapping_FeatureToMs2Indices "OpenMS::FeatureMapping::FeatureToMs2Indices":
        FeatureMapping_FeatureToMs2Indices()
        FeatureMapping_FeatureToMs2Indices(FeatureMapping_FeatureToMs2Indices)






