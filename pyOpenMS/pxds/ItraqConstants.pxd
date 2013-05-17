from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from Map cimport *
from Peak2D cimport *
from String cimport *
from StringList cimport *
from Matrix cimport *

# #ctypedef Map<Int, ChannelInfo> ChannelMapType;
# cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>" namespace "OpenMS":
#     
#     cdef cppclass ItraqConstants "OpenMS::ItraqConstants":
#         ItraqConstants(ItraqConstants) nogil except + #wrap-ignore
#         # TODO
#         # Int CHANNEL_COUNT()
#         # Int CHANNELS_FOURPLEX()
#         # Int CHANNELS_EIGHTPLEX()
#         # Int CHANNELS_TMT_SIXPLEX()
#         # double ISOTOPECORRECTIONS_FOURPLEX()
#         # double ISOTOPECORRECTIONS_EIGHTPLEX()
#         # double ISOTOPECORRECTIONS_TMT_SIXPLEX()
#         # StringList getIsotopeMatrixAsStringList(int itraq_type, IsotopeMatrices & isotope_corrections) nogil except +
#         # void updateIsotopeMatrixFromStringList(int itraq_type, StringList & channels, IsotopeMatrices & isotope_corrections) nogil except +
#         # void initChannelMap(int itraq_type, ChannelMapType & map_) nogil except +
#         # void updateChannelMap(StringList & active_channels, ChannelMapType & map_) nogil except +
#         # Matrix[ double ] translateIsotopeMatrix(int & itraq_type, IsotopeMatrices & isotope_corrections) nogil except +

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>" namespace "OpenMS::ItraqConstants":

    cdef enum ITRAQ_TYPES:

        FOURPLEX, EIGHTPLEX, TMT_SIXPLEX, SIZE_OF_ITRAQ_TYPES

    cdef cppclass ChannelInfo "OpenMS::ItraqConstants::ChannelInfo":
        ChannelInfo(ChannelInfo) nogil except + #wrap-ignore
        # TODO string variable
        # String description
        Int name
        Int id
        # TODO
        # NAMESPACE # Peak2D::CoordinateType center
        bool active

