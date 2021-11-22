from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from StringList cimport *
from AASequence cimport *
from NASequence cimport *

cdef extern from "<OpenMS/METADATA/ID/IdentifiedSequence.h>" namespace "OpenMS::IdentificationDataInternal":
    
    cdef cppclass IdentifiedSequence[SeqType]:
        # wrap-instances:
        #   IdentifiedPeptide := IdentifiedSequence[AASequence]
        #   IdentifiedOligo := IdentifiedSequence[NASequence]
        IdentifiedSequence() nogil except + # wrap-ignore
        IdentifiedSequence(SeqType& seq) nogil except +
        IdentifiedSequence(IdentifiedSequence[SeqType]&) nogil except +

        # IdentifiedSequence operator+=(const IdentifiedSequence& other) nogil except +
        bool allParentsAreDecoys() nogil except +

    cdef cppclass IdentifiedPeptideRef "OpenMS::IdentificationDataInternal::IdentifiedPeptideRef":
        # IdentifiedPeptideRef() nogil except + # wrap-ignore
        IdentifiedPeptideRef(IdentifiedPeptideRef&) nogil except + # compiler

        AASequence getAASequence() nogil except +
        IdentifiedSequence[AASequence] getIdentifiedPeptide() nogil except +

        bool operator==(IdentifiedPeptideRef) nogil except +
        # bool operator<IdentifiedPeptideRef) nogil except +
        # bool operator!=(IdentifiedPeptideRef) nogil except +

