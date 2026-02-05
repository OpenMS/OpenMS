from libcpp.vector cimport vector as libcpp_vector
from AASequence cimport AASequence

# Bindings for the OpenMS SequenceCoverage utility.
# This exposes sequence coverage calculation for proteins and peptides.
cdef extern from "<OpenMS/CHEMISTRY/SequenceCoverage.h>" namespace "OpenMS":

    cdef cppclass SequenceCoverage:
        # wrap-doc:
        #  Compute sequence coverage of a protein by peptide sequences
        # Default constructor
        SequenceCoverage() except + nogil

        # Copy constructor
        SequenceCoverage(SequenceCoverage&) except + nogil

        # Compute the percentage sequence coverage of a protein
        # by a given list of peptide sequences.
        @staticmethod
        double getCoverage(
            const AASequence& protein,
            const libcpp_vector[AASequence]& peptides
        ) except + nogil
