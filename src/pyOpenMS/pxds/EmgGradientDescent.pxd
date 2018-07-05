from MSChromatogram cimport *
from MSSpectrum cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/MATH/MISC/EmgGradientDescent.h>" namespace "OpenMS":

    cdef cppclass EmgGradientDescent(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        EmgGradientDescent() nogil except +
        EmgGradientDescent(EmgGradientDescent) nogil except +

        void getDefaultParameters(Param&) nogil except +

        void fitEMGPeakModel(MSChromatogram& input_peak, MSChromatogram& output_peak) nogil except +
        void fitEMGPeakModel(MSSpectrum& input_peak, MSSpectrum& output_peak) nogil except +
        void fitEMGPeakModel(MSChromatogram& input_peak, MSChromatogram& output_peak, double left_pos, double right_pos) nogil except +
        void fitEMGPeakModel(MSSpectrum& input_peak, MSSpectrum& output_peak, double left_pos, double right_pos) nogil except +
