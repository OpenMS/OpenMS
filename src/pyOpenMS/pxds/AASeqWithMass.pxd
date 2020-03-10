from Types cimport *
from OPXLDataStructs cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS":

    cdef cppclass AASeqWithMass "OpenMS::OPXLDataStructs::AASeqWithMass":

        AASeqWithMass(AASeqWithMass) nogil except +
        AASeqWithMass() nogil except +

        double peptide_mass
        AASequence peptide_seq
        PeptidePosition position
