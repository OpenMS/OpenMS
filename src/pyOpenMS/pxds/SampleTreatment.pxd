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

        SampleTreatment(SampleTreatment &) except + nogil  # wrap-doc:Base class for sample treatments (Digestion, Modification, Tagging, ...)
        SampleTreatment(const String & type_) except + nogil 

        bool operator==(SampleTreatment & rhs) except + nogil 
        String getType() except + nogil  # wrap-doc:Returns the treatment type
        String getComment() except + nogil  # wrap-doc:Returns the description of the sample treatment
        void setComment(const String & comment) except + nogil  # wrap-doc:Sets the description of the sample treatment
        # POINTER # SampleTreatment * clone() except + nogil 
