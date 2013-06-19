from Types cimport *
from Feature cimport *
from DocumentIdentifier cimport *
from ProteinIdentification cimport *
from DataProcessing cimport *

cdef extern from "<OpenMS/KERNEL/FeatureMap.h>" namespace "OpenMS":
    
    cdef cppclass AnnotationStatistics "OpenMS::AnnotationStatistics":
        AnnotationStatistics() nogil except +
        AnnotationStatistics(AnnotationStatistics) nogil except +
        libcpp_vector[ size_t ] states
        bool operator==(AnnotationStatistics & rhs) nogil except +
        # NAMESPACE # AnnotationStatistics  operator+=(BaseFeature::AnnotationState state) nogil except +

