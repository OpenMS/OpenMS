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
        # no-pxd-import

        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricQuantitationMethod() except + nogil 
        IsobaricQuantitationMethod(IsobaricQuantitationMethod &) except + nogil 
        String getName() except + nogil 
        libcpp_vector[IsobaricChannelInformation]  getChannelInformation() except + nogil 
        Size getNumberOfChannels() except + nogil 
        Matrix[ double ] getIsotopeCorrectionMatrix() except + nogil 
        Size getReferenceChannel() except + nogil 

