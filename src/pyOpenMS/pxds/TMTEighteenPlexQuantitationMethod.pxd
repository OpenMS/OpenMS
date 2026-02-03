from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTEighteenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-doc:
        #  TMT 18-plex isobaric labeling quantitation method. Provides channel information
        #  and isotope correction matrix for Thermo Scientific TMTpro 18-plex reagents
        #  which extends TMTpro 16-plex with two additional channels

        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTEighteenPlexQuantitationMethod() except + nogil  # wrap-doc:Default constructor
        TMTEighteenPlexQuantitationMethod(TMTEighteenPlexQuantitationMethod &) except + nogil  # wrap-doc:Copy constructor
