from Types cimport *
from Param cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationStandards.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationStandards:

        AbsoluteQuantitationStandards() nogil except +
        AbsoluteQuantitationStandards(AbsoluteQuantitationStandards)  nogil except + #wrap-ignore
