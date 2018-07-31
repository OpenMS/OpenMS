from Types cimport *
from FeatureMap cimport *
from MRMFeatureQC cimport *
from DefaultParamHandler cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureFilter(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        MRMFeatureFilter() nogil except +
        MRMFeatureFilter(MRMFeatureFilter &) nogil except +

        void FilterFeatureMap(FeatureMap features, MRMFeatureQC filter_criteria, TargetedExperiment transitions) nogil except +

        # could add support for other members if needed
        # however, they are only used internally for now
