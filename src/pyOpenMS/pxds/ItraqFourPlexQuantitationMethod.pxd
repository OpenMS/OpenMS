from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass ItraqFourPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        ItraqFourPlexQuantitationMethod() nogil except +
        ItraqFourPlexQuantitationMethod(ItraqFourPlexQuantitationMethod) nogil except +

