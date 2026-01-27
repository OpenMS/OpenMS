from libcpp.vector cimport vector
from pyopenms cimport AASequence

cdef extern from "<OpenMS/CHEMISTRY/SequenceCoverage.h>" namespace "OpenMS":

    cdef cppclass SequenceCoverage:
        @staticmethod
        double getCoverage(
            const AASequence& protein,
            const vector[AASequence]& peptides
        ) except +
