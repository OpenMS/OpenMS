from Types cimport *
from DefaultParamHandler cimport *
from Matrix cimport *
from String cimport *
from Peak2D cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>" namespace "OpenMS::IsobaricQuantitationMethod":
    
    cdef cppclass IsobaricChannelInformation "OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation":
        IsobaricChannelInformation(IsobaricChannelInformation) nogil except + #wrap-ignore
        Int name
        Int id
        String description
        double center
        IsobaricChannelInformation(Int name, Int id_, String description, double center) nogil except +

