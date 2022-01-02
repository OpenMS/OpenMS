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

        SampleTreatment(SampleTreatment &) nogil except + # wrap-doc:Base class for sample treatments (Digestion, Modification, Tagging, ...)
        SampleTreatment(const String & type_) nogil except +

        bool operator==(SampleTreatment & rhs) nogil except +
        String getType() nogil except + # wrap-doc:Returns the treatment type
        String getComment() nogil except + # wrap-doc:Returns the description of the sample treatment
        void setComment(const String & comment) nogil except + # wrap-doc:Sets the description of the sample treatment
        # POINTER # SampleTreatment * clone() nogil except +
