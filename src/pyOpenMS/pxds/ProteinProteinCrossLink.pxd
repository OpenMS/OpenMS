from Types cimport *
from libcpp cimport bool
from PeptideHit cimport *
from MSExperiment cimport *
from ResidueModification cimport *
from FASTAFile cimport *
from ProteaseDigestion cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
    
    cdef cppclass ProteinProteinCrossLink "OpenMS::OPXLDataStructs::ProteinProteinCrossLink":
        ProteinProteinCrossLink(ProteinProteinCrossLink) nogil except + #wrap-ignore
        AASequence alpha
        AASequence beta
        libcpp_pair[ ptrdiff_t, ptrdiff_t] cross_link_position
        double cross_linker_mass
        String cross_linker_name
        TermSpecificity term_spec_alpha
        TermSpecificity term_spec_beta
        # ProteinProteinCrossLinkType getType() nogil except +
        bool operator==(ProteinProteinCrossLink & other) nogil except +

