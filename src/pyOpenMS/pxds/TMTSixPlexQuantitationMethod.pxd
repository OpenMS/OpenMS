from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTSixPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-doc:
        #  TMT 6-plex isobaric labeling quantitation method. Provides channel information
        #  and isotope correction matrix for Thermo Scientific TMT 6-plex reagents
        #  with reporter ions at 126, 127, 128, 129, 130, and 131 m/z

        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTSixPlexQuantitationMethod() except + nogil  # wrap-doc:Default constructor
        TMTSixPlexQuantitationMethod(TMTSixPlexQuantitationMethod &) except + nogil  # wrap-doc:Copy constructor
