from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass TMTSixteenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTSixteenPlexQuantitationMethod() nogil except +
        TMTSixteenPlexQuantitationMethod(TMTSixteenPlexQuantitationMethod) nogil except +
