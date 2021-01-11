from Types cimport *
from libcpp cimport bool
from String cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/SampleTreatment.h>" namespace "OpenMS":
    
    cdef cppclass SampleTreatment(MetaInfoInterface) :
        # wrap-ignore
        # no-pxd-import
        # ABSTRACT class
        #
        # wrap-inherits:
        #  MetaInfoInterface

        SampleTreatment(SampleTreatment) nogil except +
        SampleTreatment(const String & type_) nogil except +

        bool operator==(SampleTreatment & rhs) nogil except +
        String getType() nogil except +
        String getComment() nogil except +
        void setComment(const String & comment) nogil except +
        # POINTER # SampleTreatment * clone() nogil except +


