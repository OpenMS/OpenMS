from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass ItraqFourPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        ItraqFourPlexQuantitationMethod() nogil except + # wrap-doc:iTRAQ 4 plex quantitation to be used with the IsobaricQuantitation
        ItraqFourPlexQuantitationMethod(ItraqFourPlexQuantitationMethod &) nogil except +


