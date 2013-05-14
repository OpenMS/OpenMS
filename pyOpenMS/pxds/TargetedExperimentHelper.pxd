from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair
from libcpp cimport bool as libcpp_pair

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass RetentionTime:
        RetentionTime() nogil except +
        RetentionTime(RetentionTime) nogil except +
        String software_ref

    cdef cppclass Protein:

        Protein() nogil except +
        Protein(Protein) nogil except +

    cdef cppclass Peptide: # TODO add CVTermList

        Peptide() nogil except +
        Peptide(Peptide) nogil except +
        double getRetentionTime() nogil except +
        String sequence
        libcpp_vector[RetentionTime] rts

