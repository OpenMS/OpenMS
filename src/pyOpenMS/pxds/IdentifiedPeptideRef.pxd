from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from StringList cimport *

cdef extern from "<OpenMS/METADATA/ID/IdentifiedSequence.h>" namespace "OpenMS::IdentificationDataInternal":
    
    cdef cppclass IdentifiedPeptideRef "OpenMS::IdentificationDataInternal::IdentifiedPeptideRef":
        IdentifiedPeptideRef() nogil except + # wrap-ignore
        IdentifiedPeptideRef(IdentifiedPeptideRef&) nogil except + # compiler
