from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair
from libcpp cimport bool as libcpp_pair

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass Protein:

        Protein() nogil except +
        Protein(Protein) nogil except +

    cdef cppclass Peptide:

        Peptide() nogil except +
        Peptide(Peptide) nogil except +
        double getRetentionTime() except +

