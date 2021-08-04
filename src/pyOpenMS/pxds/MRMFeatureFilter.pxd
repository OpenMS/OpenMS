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
        MRMFeatureFilter(MRMFeatureFilter &) nogil except + # compiler

        void FilterFeatureMap(FeatureMap features, MRMFeatureQC filter_criteria, TargetedExperiment transitions) nogil except +
            # wrap-doc:
                #   Flags or filters features and subordinates in a FeatureMap
                #   -----
                #   :param features: FeatureMap to flag or filter
                #   :param filter_criteria: MRMFeatureQC class defining QC parameters
                #   :param transitions: Transitions from a TargetedExperiment

        # could add support for other members if needed
        # however, they are only used internally for now
