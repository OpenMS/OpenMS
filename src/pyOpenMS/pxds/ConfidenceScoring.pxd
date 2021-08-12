from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h>" namespace "OpenMS":

    cdef cppclass ConfidenceScoring:

        ConfidenceScoring() nogil except +
        ConfidenceScoring(ConfidenceScoring &) nogil except + # compiler

        void initialize(TargetedExperiment & targeted, Size n_decoys, Size n_transitions, TransformationDescription trafo) nogil except +
        void initializeGlm(double intercept, double rt_coef, double int_coef) nogil except +
        void scoreMap(FeatureMap & map) nogil except + # wrap-doc:Score a feature map -> make sure the class is properly initialized

