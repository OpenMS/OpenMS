from Types cimport *
from Feature cimport *
from DocumentIdentifier cimport *
from ProteinIdentification cimport *
from DataProcessing cimport *

cdef extern from "<OpenMS/KERNEL/FeatureMap.h>" namespace "OpenMS":
    
    cdef cppclass AnnotationStatistics "OpenMS::AnnotationStatistics":
        AnnotationStatistics() except + nogil 
        AnnotationStatistics(AnnotationStatistics &) except + nogil 
        libcpp_vector[ size_t ] states
        bool operator==(AnnotationStatistics & rhs) except + nogil 
        # NAMESPACE # AnnotationStatistics  operator+=(BaseFeature::AnnotationState state) except + nogil 

