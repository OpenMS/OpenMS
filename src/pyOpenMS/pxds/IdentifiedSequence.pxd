from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
from ScoredProcessingResult cimport *
from ParentMatch cimport *
from AASequence cimport *
from NASequence cimport *


cdef extern from "<OpenMS/METADATA/ID/IdentifiedSequence.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass IdentifiedSequence[SeqType](ScoredProcessingResult):

    IdentifiedSequence() nogil except +

    IdentifiedSequence(IdentifiedSequence & other) nogil except +

    SeqType sequence

    ParentMatches parent_matches

    IdentifiedSequence(const SeqType & sequence, const ParentMatches & parent_matchs, const AppliedProcessingSteps & steps_and_scores) # Note that default args don't work

    IdentifiedSequence & merge(const IdentifiedSequence & other) nogil except +

    bool allParentsAreDecoys() nogil except +

  #ctypedef IdentifiedSequence[AASequence] IdentifiedPeptide # WHY DOES THIS NEED A CONVERTER FIXME 
  #ctypedef IdentifiedSequence[NASequence] IdentifiedOligo

  cdef cppclass IdentifiedPeptides:
    IdentifiedPeptides() nogil except +
    IdentifiedPeptides(IdentifiedPeptides & other) nogil except +
    #IdentifiedSequence operator[](size_t index) #wrap-upper-limit:size() #TODO: Add some sort of access to get the IdentifiedPeptides back out

  cdef cppclass IdentifiedOligos:
    IdentifiedOligos() nogil except +
    IdentifiedOligos(IdentifiedOligos & other) nogil except +
    
  cdef cppclass IdentifiedPeptideRef:
    IdentifiedPeptideRef() nogil except +
    IdentifiedPeptideRef(const IdentifiedPeptideRef & other) nogil except +
    bool operator!=(const IdentifiedPeptideRef & other) nogil except +
    bool operator<(const IdentifiedPeptideRef & other) nogil except +
    #IdentifiedPeptide deref()

  cdef cppclass IdentifiedOligoRef:
    IdentifiedOligoRef() nogil except +
    IdentifiedOligoRef(const IdentifiedOligoRef & other) nogil except +
    bool operator!=(const IdentifiedOligoRef & other) nogil except +
    bool operator<(const IdentifiedOligoRef & other) nogil except +
    #IdentifiedOligo deref()