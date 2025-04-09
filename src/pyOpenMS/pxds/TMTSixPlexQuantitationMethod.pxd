from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass TMTSixPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTSixPlexQuantitationMethod() except + nogil 
        TMTSixPlexQuantitationMethod(TMTSixPlexQuantitationMethod &) except + nogil 
