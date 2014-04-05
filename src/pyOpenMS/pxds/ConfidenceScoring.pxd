from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from FeatureMap cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h>" namespace "OpenMS":

    cdef cppclass ConfidenceScoring:

        ConfidenceScoring() nogil except +
        ConfidenceScoring(ConfidenceScoring) nogil except +

        void initialize(TargetedExperiment & targeted, Size n_decoys, Size n_transitions, TransformationDescription trafo) nogil except +
        void initializeGlm(double intercept, double rt_coef, double int_coef) nogil except +
        void scoreMap(FeatureMap[Feature] & map) nogil except +

