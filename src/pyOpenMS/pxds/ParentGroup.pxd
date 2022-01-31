from MetaData cimport *
from MetaInfoInterface cimport *
from ScoredProcessingResult cimport *
from String cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
from libcpp.vector cimport vector as libcpp_vector
from ParentSequence cimport *
from ScoreType cimport *

cdef extern from "<OpenMS/METADATA/ID/ParentGroup.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass ParentGroup(MetaInfoInterface):
    ParentGroup() nogil except +
    libcpp_map[IteratorWrapper[ setSTit, ScoreType ], double] scores
    libcpp_set[ParentSequenceRef] parent_refs
  
  cdef cppclass ParentGroups:
    ParentGroups() nogil except +
    ParentGroups(ParentGroups other) nogil except +

  cdef cppclass ParentGroupRef:
    ParentGroupRef() nogil except +
    ParentGroupRef(const ParentGroupRef  other) nogil except +
    bool operator!=(const ParentGroupRef  other) nogil except +
    bool operator<(const ParentGroupRef  other) nogil except +
    ParentSequence deref()

  cdef cppclass ParentGroupSet(ScoredProcessingResult):
    ParentGroupSet() nogil except +
    String label
    ParentGroups groups
    
  ctypedef libcpp_vector[ParentGroupSet] ParentGroupSets