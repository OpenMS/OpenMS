from libcpp.vector cimport vector
from AASequence cimport AASequence

cdef extern from "<OpenMS/DATASTRUCTURES/ListUtils.h>" namespace "OpenMS":
    ctypedef vector[AASequence] AASequenceList