from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTElevenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTElevenPlexQuantitationMethod() nogil except +
        TMTElevenPlexQuantitationMethod(TMTElevenPlexQuantitationMethod) nogil except +

