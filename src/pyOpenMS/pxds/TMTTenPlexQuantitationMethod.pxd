from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTTenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-doc:
        #  TMT 10-plex isobaric labeling quantitation method. Provides channel information
        #  and isotope correction matrix for Thermo Scientific TMT 10-plex reagents
        #  with reporter ions spanning the 126-131 m/z range with additional mass variants

        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTTenPlexQuantitationMethod() except + nogil  # wrap-doc:Default constructor
        TMTTenPlexQuantitationMethod(TMTTenPlexQuantitationMethod &) except + nogil  # wrap-doc:Copy constructor
