from Types cimport *
from OPXLDataStructs cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS":

    cdef cppclass AASeqWithMass "OpenMS::OPXLDataStructs::AASeqWithMass":

        AASeqWithMass() nogil except +
        AASeqWithMass(AASeqWithMass &) nogil except + # compiler

        double peptide_mass
        AASequence peptide_seq
        PeptidePosition position
