from Types cimport *
from DefaultParamHandler cimport *
from Matrix cimport *
from String cimport *
from Peak2D cimport *
from IsobaricChannelInformation cimport *

# typedef std::vector<IsobaricChannelInformation> IsobaricChannelList;
cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricQuantitationMethod(DefaultParamHandler) :
        # wrap-ignore
        # ABSTRACT class

        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricQuantitationMethod() nogil except +
        IsobaricQuantitationMethod(IsobaricQuantitationMethod) nogil except + #wrap-ignore
        # String  getName() nogil except +
        libcpp_vector[IsobaricChannelInformation]  getChannelInformation() nogil except +
        Size getNumberOfChannels() nogil except +
        Matrix[ double ] getIsotopeCorrectionMatrix() nogil except +
        Size getReferenceChannel() nogil except +

