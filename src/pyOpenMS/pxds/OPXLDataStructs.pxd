from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS":
    
    cdef cppclass OPXLDataStructs "OpenMS::OPXLDataStructs":
        OPXLDataStructs(OPXLDataStructs) nogil except + #wrap-ignore

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
    cdef enum ProteinProteinCrossLinkType "OpenMS::OPXLDataStructs::ProteinProteinCrossLinkType":
        #wrap-attach:
        #    OPXLDataStructs
        CROSS
        MONO
        LOOP
        NUMBER_OF_CROSS_LINK_TYPES

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
    cdef enum PeptidePosition "OpenMS::OPXLDataStructs::PeptidePosition":
        #wrap-attach:
        #    OPXLDataStructs
        INTERNAL
        C_TERM
        N_TERM

