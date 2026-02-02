from AASequence cimport AASequence
from AASequenceList cimport AASequenceList

# Bindings for the OpenMS SequenceCoverage utility.
# This exposes sequence coverage calculation for proteins and peptides.
cdef extern from "<OpenMS/CHEMISTRY/SequenceCoverage.h>" namespace "OpenMS":

    cdef cppclass SequenceCoverage:
        # Default constructor
        SequenceCoverage() except + nogil

        # Copy constructor
        SequenceCoverage(SequenceCoverage&) except + nogil

        # Compute the percentage sequence coverage of a protein
        # by a given list of peptide sequences.
        @staticmethod
        double getCoverage(
            const AASequence& protein,
            const AASequenceList& peptides
        ) except + nogil
