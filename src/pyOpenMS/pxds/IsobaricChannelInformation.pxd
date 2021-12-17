from Types cimport *
from DefaultParamHandler cimport *
from Matrix cimport *
from String cimport *
from Peak2D cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>" namespace "OpenMS::IsobaricQuantitationMethod":

    cdef cppclass IsobaricChannelInformation "OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation":

        IsobaricChannelInformation(String name, Int id_, String description, double center, Int channel_id_minus_2, Int channel_id_minus_1, Int channel_id_plus_1, Int channel_id_plus_2) nogil except +
        IsobaricChannelInformation(IsobaricChannelInformation &) nogil except +
        String name
        Int id
        String description
        double center
        Int channel_id_minus_2
        Int channel_id_minus_1
        Int channel_id_plus_1
        Int channel_id_plus_2
