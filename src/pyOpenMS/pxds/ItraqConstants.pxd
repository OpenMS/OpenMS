from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from Map cimport *
from Peak2D cimport *
from String cimport *
from StringList cimport *
from Matrix cimport *

# ctypedef Map<Int, ChannelInfo> ChannelMapType; 
# ctypedef std::vector<Matrix<double> > IsotopeMatrices;
cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>" namespace "OpenMS":
    
    cdef cppclass ItraqConstants "OpenMS::ItraqConstants":
        # wrap-doc:
                #   Some constants used throughout iTRAQ classes
                #   -----
                #   Constants for iTRAQ experiments and a ChannelInfo structure to store information about a single channel

        ItraqConstants() nogil except + # compiler
        ItraqConstants(ItraqConstants &) nogil except + # compiler

        # Int CHANNEL_COUNT()
        # Int CHANNELS_FOURPLEX()
        # Int CHANNELS_EIGHTPLEX()
        # Int CHANNELS_TMT_SIXPLEX()
        # double ISOTOPECORRECTIONS_FOURPLEX()
        # double ISOTOPECORRECTIONS_EIGHTPLEX()
        # double ISOTOPECORRECTIONS_TMT_SIXPLEX()

        StringList getIsotopeMatrixAsStringList(int itraq_type, libcpp_vector[Matrix[double] ] & isotope_corrections) nogil except +
            # wrap-doc:
                #   Convert isotope correction matrix to stringlist
                #   -----
                #   Each line is converted into a string of the format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
                #   Useful for creating parameters or debug output
                #   -----
                #   :param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
                #   :param isotope_corrections: Vector of the two matrices (4plex, 8plex)

        void updateIsotopeMatrixFromStringList(int itraq_type, StringList & channels, libcpp_vector[Matrix[double] ] & isotope_corrections) nogil except +
            # wrap-doc:
                #   Convert strings to isotope correction matrix rows
                #   -----
                #   Each string of format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
                #   is parsed and the corresponding channel(row) in the matrix is updated
                #   Not all channels need to be present, missing channels will be left untouched
                #   Useful to update the matrix with user isotope correction values
                #   -----
                #   :param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
                #   :param channels: New channel isotope values as strings
                #   :param isotope_corrections: Vector of the two matrices (4plex, 8plex)

        # void initChannelMap(int itraq_type, ChannelMapType & map_) nogil except +
        # void updateChannelMap(StringList & active_channels, ChannelMapType & map_) nogil except +
        Matrix[ double ] translateIsotopeMatrix(int & itraq_type, libcpp_vector[Matrix[double] ] & isotope_corrections) nogil except +

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>" namespace "OpenMS::ItraqConstants":

    cdef enum ITRAQ_TYPES:

        FOURPLEX, EIGHTPLEX, TMT_SIXPLEX, SIZE_OF_ITRAQ_TYPES

    cdef cppclass ChannelInfo "OpenMS::ItraqConstants::ChannelInfo":
        ChannelInfo(ChannelInfo) nogil except + #wrap-ignore
        libcpp_string description
        Int name
        Int id
        double center
        bool active

