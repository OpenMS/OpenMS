from Types cimport *
# from IsobaricQuantitationMethod cimport *
from DefaultParamHandler cimport *
from Matrix cimport *
from String cimport *
from Peak2D cimport *
from IsobaricChannelInformation cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTThirtyTwoPlexQuantitationMethod(DefaultParamHandler) :
        # wrap-inherits:
        # IsobaricQuantitationMethod
        TMTThirtyTwoPlexQuantitationMethod() except + nogil
        TMTThirtyTwoPlexQuantitationMethod(TMTThirtyTwoPlexQuantitationMethod &) except + nogil
        String getMethodName() except + nogil 
        libcpp_vector[IsobaricChannelInformation]  getChannelInformation() except + nogil 
        Size getNumberOfChannels() except + nogil 
        Matrix[ double ] getIsotopeCorrectionMatrix() except + nogil 
        Size getReferenceChannel() except + nogil 
