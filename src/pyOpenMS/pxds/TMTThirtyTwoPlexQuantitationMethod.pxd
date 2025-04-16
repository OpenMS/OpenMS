from Types cimport *
# from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTThirtyTwoPlexQuantitationMethod(DefaultParamHandler) :
        # wrap-inherits:
        #   IsobaricQuantitationMethod
        TMTThirtyTwoPlexQuantitationMethod() except + nogil
        TMTThirtyTwoPlexQuantitationMethod(TMTThirtyTwoPlexQuantitationMethod &) except + nogil
