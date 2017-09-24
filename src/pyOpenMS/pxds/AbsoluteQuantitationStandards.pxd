from Types cimport *
from Param cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationStandards.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationStandards:
        # wrap-ignore

        AbsoluteQuantitationStandards() nogil except +