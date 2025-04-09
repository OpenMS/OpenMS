from Types cimport *
from BaseFeature cimport *
from MSExperiment cimport *
from KDTreeFeatureMaps cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>" namespace "OpenMS":

    cdef cppclass FeatureMapping "OpenMS::FeatureMapping":
       FeatureMapping() except + nogil  # compiler
       FeatureMapping(FeatureMapping &) except + nogil  # compiler

    cdef cppclass FeatureMapping_FeatureMappingInfo "OpenMS::FeatureMapping::FeatureMappingInfo":
       FeatureMapping_FeatureMappingInfo() except + nogil 
       FeatureMapping_FeatureMappingInfo(FeatureMapping_FeatureMappingInfo &) except + nogil 

    cdef cppclass FeatureMapping_FeatureToMs2Indices "OpenMS::FeatureMapping::FeatureToMs2Indices":
        FeatureMapping_FeatureToMs2Indices() except + nogil 
        FeatureMapping_FeatureToMs2Indices(FeatureMapping_FeatureToMs2Indices &) except + nogil 
        
cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>" namespace "OpenMS::FeatureMapping":

       # wrap static method:
       FeatureMapping_FeatureToMs2Indices assignMS2IndexToFeature(MSExperiment& spectra,
                                                                  FeatureMapping_FeatureMappingInfo& fm_info,
                                                                  double precursor_mz_tolerance,
                                                                  double precursor_rt_tolerance,
                                                                  bool ppm) except + nogil   # wrap-attach:FeatureMapping
