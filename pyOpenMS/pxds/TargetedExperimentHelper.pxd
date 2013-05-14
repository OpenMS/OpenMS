#from DataValue cimport *
#from String cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>" namespace "OpenMS::TargetedExperimentHelper":

    cdef cppclass Protein:

        Protein() nogil except +
        Protein(Protein) nogil except +

    cdef cppclass Peptide:

        Peptide() nogil except +
        Peptide(Peptide) nogil except +

