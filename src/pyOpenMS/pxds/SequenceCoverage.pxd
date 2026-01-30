from libcpp.vector cimport vector
from pyopenms cimport AASequence

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
            const vector[AASequence]& peptides
        ) except + nogil