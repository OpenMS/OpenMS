from OPXLDataStructs cimport *
from String cimport String
from libcpp cimport bool
from ResidueModification cimport TermSpecificity
from AASequence cimport AASequence


cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":

    cdef cppclass ProteinProteinCrossLink "OpenMS::OPXLDataStructs::ProteinProteinCrossLink":

        ProteinProteinCrossLink(ProteinProteinCrossLink) nogil except +
        ProteinProteinCrossLink() nogil except +


        AASequence alpha
        AASequence beta
        libcpp_pair[ ptrdiff_t, ptrdiff_t] cross_link_position
        double cross_linker_mass
        String cross_linker_name
        TermSpecificity term_spec_alpha
        TermSpecificity term_spec_beta
        int precursor_correction

        ProteinProteinCrossLinkType getType() nogil except +
        bool operator==(ProteinProteinCrossLink other) nogil except +
