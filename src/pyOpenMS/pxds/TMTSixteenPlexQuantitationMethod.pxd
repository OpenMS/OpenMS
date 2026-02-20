from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTSixteenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-doc:
        #  TMT 16-plex isobaric labeling quantitation method. Provides channel information
        #  and isotope correction matrix for Thermo Scientific TMTpro 16-plex reagents
        #  with reporter ions spanning the 126-134 m/z range

        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTSixteenPlexQuantitationMethod() except + nogil  # wrap-doc:Default constructor
        TMTSixteenPlexQuantitationMethod(TMTSixteenPlexQuantitationMethod &) except + nogil  # wrap-doc:Copy constructor
