from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTThirtyFivePlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTThirtyFivePlexQuantitationMethod(IsobaricQuantitationMethod) :
        # IsobaricQuantitationMethod
        TMTThirtyFivePlexQuantitationMethod() except + nogil
        TMTThirtyFivePlexQuantitationMethod(TMTThirtyFivePlexQuantitationMethod &) except + nogil