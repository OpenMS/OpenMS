from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTThirtyFivePlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTThirtyFivePlexQuantitationMethod(IsobaricQuantitationMethod) :
        # IsobaricQuantitationMethod
        TMTThirtyFIvePlexQuantitationMethod() except + nogil
        TMTThirtyFIvePlexQuantitationMethod(TMTThirtyFivePlexQuantitationMethod &) except + nogil