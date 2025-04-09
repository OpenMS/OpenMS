from Types cimport *
from DefaultParamHandler cimport *
from Matrix cimport *
from String cimport *
from Peak2D cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>" namespace "OpenMS::IsobaricQuantitationMethod":

    cdef cppclass IsobaricChannelInformation "OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation":

        IsobaricChannelInformation(String name, Int id_, String description, double center, libcpp_vector[int] affected_channels) except + nogil 
        IsobaricChannelInformation(IsobaricChannelInformation &) except + nogil 
        String name
        Int id
        String description
        double center
        libcpp_vector[int] affected_channels
