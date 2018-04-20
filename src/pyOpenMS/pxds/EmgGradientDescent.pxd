from MSChromatogram cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/MATH/MISC/EmgGradientDescent.h>" namespace "OpenMS":

    cdef cppclass EmgGradientDescent(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        EmgGradientDescent() nogil except +
        EmgGradientDescent(EmgGradientDescent) nogil except +

        void getDefaultParameters(Param&) nogil except +

        void fitEMGPeakModel(const MSChromatogram& input_peak, MSChromatogram& output_peak) nogil except +
        void fitEMGPeakModel(const MSSpectrum& input_peak, MSSpectrum& output_peak) nogil except +
