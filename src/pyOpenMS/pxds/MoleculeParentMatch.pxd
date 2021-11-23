from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/METADATA/ID/MoleculeParentMatch.h>" namespace "OpenMS::IdentificationDataInternal":
    
    cdef cppclass MoleculeParentMatch "OpenMS::IdentificationDataInternal::MoleculeParentMatch":
        MoleculeParentMatch() nogil except +
        MoleculeParentMatch(MoleculeParentMatch&) nogil except +

        MoleculeParentMatch(Size start_pos,
                            Size end_pos,
                            String left_neighbor,
                            String right_neighbor) nogil except +

        # members
        Size start_pos
        Size end_pos

        String left_neighbor
        String right_neighbor
