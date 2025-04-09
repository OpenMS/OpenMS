from MSChromatogram cimport *
from MSSpectrum cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/MATH/MISC/EmgGradientDescent.h>" namespace "OpenMS":

    cdef cppclass EmgGradientDescent(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        EmgGradientDescent() except + nogil  # wrap-doc:Compute the area, background and shape metrics of a peak
        EmgGradientDescent(EmgGradientDescent &) except + nogil  # compiler

        void getDefaultParameters(Param&) except + nogil 

        void fitEMGPeakModel(MSChromatogram& input_peak, MSChromatogram& output_peak) except + nogil 
        void fitEMGPeakModel(MSSpectrum& input_peak, MSSpectrum& output_peak) except + nogil 
        void fitEMGPeakModel(MSChromatogram& input_peak, MSChromatogram& output_peak, double left_pos, double right_pos) except + nogil 
        void fitEMGPeakModel(MSSpectrum& input_peak, MSSpectrum& output_peak, double left_pos, double right_pos) except + nogil 
