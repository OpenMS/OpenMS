from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from FeatureMap cimport *
from MRMFeatureQC import *
from TargetedExperiment import *

cdef extern from "<OpenMS/FORMAT/MRMFeatureFilter.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureFilter:

        MRMFeatureFilter() nogil except +
        MRMFeatureFilter(MRMFeatureFilter &) nogil except +

        void FilterFeatureMap(FeatureMap features, MRMFeatureQC filter_criteria, TargetedExperiment transitions) nogil except +

        # could add support for other members if needed
        # however, they are only used internally for now