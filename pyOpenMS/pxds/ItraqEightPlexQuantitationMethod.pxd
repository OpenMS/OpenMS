from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass ItraqEightPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        ItraqEightPlexQuantitationMethod() nogil except +
        ItraqEightPlexQuantitationMethod(ItraqEightPlexQuantitationMethod) nogil except +

