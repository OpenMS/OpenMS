from libcpp.vector cimport vector
from pyopenms cimport AASequence

# Bindings for the OpenMS SequenceCoverage utility.
# Exposes C++ sequence coverage computation to pyOpenMS.
cdef extern from "<OpenMS/CHEMISTRY/SequenceCoverage.h>" namespace "OpenMS":

    cdef cppclass SequenceCoverage:
        SequenceCoverage() except + nogil
        SequenceCoverage(SequenceCoverage&) except + nogil
        # Compute the percentage sequence coverage of a protein by given peptides.
        @staticmethod
        double getCoverage(
            const AASequence& protein,
            const vector[AASequence]& peptides
        ) except +
