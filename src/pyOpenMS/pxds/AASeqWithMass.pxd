from OPXLDataStructs cimport *
from AASequence cimport AASequence

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS":

    cdef cppclass AASeqWithMass "OpenMS::OPXLDataStructs::AASeqWithMass":
        AASeqWithMass(AASeqWithMass) nogil except +
        #wrap-attach:
        #    OPXLDataStructs

        double peptide_mass
        AASequence peptide_seq
        PeptidePosition position
