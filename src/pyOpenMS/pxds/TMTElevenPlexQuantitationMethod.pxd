from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTElevenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-doc:
        #  TMT 11-plex isobaric labeling quantitation method. Provides channel information
        #  and isotope correction matrix for Thermo Scientific TMT 11-plex reagents
        #  which extends TMT 10-plex with an additional 131C channel

        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTElevenPlexQuantitationMethod() except + nogil  # wrap-doc:Default constructor
        TMTElevenPlexQuantitationMethod(TMTElevenPlexQuantitationMethod &) except + nogil  # wrap-doc:Copy constructor
