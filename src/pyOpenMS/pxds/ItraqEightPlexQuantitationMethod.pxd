from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass ItraqEightPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        ItraqEightPlexQuantitationMethod() nogil except + # wrap-doc:iTRAQ 8 plex quantitation to be used with the IsobaricQuantitation
        ItraqEightPlexQuantitationMethod(ItraqEightPlexQuantitationMethod &) nogil except +


