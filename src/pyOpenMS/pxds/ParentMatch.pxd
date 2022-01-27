from MetaData cimport *
from MetaInfoInterface cimport *
from ScoredProcessingResult cimport *
from String cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.set cimport set as libcpp_set
from ParentSequence cimport *

cdef extern from "<OpenMS/METADATA/ID/ParentMatch.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass ParentMatch(MetaInfoInterface):
    Size start_pos
    Size end_pos

    String left_neighbor
    String right_neighbor

    #Size UNKNOWN_POSITION = Size(-1)
    #char UNKNOWN_NEIGHBOR = 'X'
    #char LEFT_TERMINUS = '['
    #char RIGHT_TERMINUS = ']'

    ParentMatch(Size start_pos, Size end_pos, String left_neighbor, String right_neighbor) #FIXME can we add any of the defaults here

    bool operator<(const ParentMatch & other)
    bool operator==(const ParentMatch & other)
    bool hasValidPositions(Size molecule_length, Size parent_length)
    
  ctypedef libcpp_map[ParentSequenceRef, libcpp_set[ParentMatch]] ParentMatches