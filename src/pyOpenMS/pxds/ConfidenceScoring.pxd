from libcpp.vector cimport vector as libcpp_vector
from FeatureMap cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h>" namespace "OpenMS":

    cdef cppclass ConfidenceScoring:

        ConfidenceScoring() except + nogil 
        ConfidenceScoring(ConfidenceScoring &) except + nogil  # compiler

        void initialize(TargetedExperiment & targeted, Size n_decoys, Size n_transitions, TransformationDescription trafo) except + nogil 
        void initializeGlm(double intercept, double rt_coef, double int_coef) except + nogil 
        void scoreMap(FeatureMap & map) except + nogil  # wrap-doc:Score a feature map -> make sure the class is properly initialized

