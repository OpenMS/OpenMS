from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>" namespace "OpenMS":
    
    cdef cppclass ItraqEightPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        ItraqEightPlexQuantitationMethod() nogil except +
        ItraqEightPlexQuantitationMethod(ItraqEightPlexQuantitationMethod) nogil except +
        # String  getName() nogil except +
        # IsobaricChannelList  getChannelInformation() nogil except +
        # Size getNumberOfChannels() nogil except +
        # Matrix[ double ] getIsotopeCorrectionMatrix() nogil except +
        # Size getReferenceChannel() nogil except +

