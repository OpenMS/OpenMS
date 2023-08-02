from Types cimport *
# from MathFunctions cimport *
from ContinuousWaveletTransform cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>" namespace "OpenMS":
    
    cdef cppclass ContinuousWaveletTransformNumIntegration(ContinuousWaveletTransform) :
        # wrap-inherits:
        #  ContinuousWaveletTransform
        ContinuousWaveletTransformNumIntegration() except + nogil 
        ContinuousWaveletTransformNumIntegration(ContinuousWaveletTransformNumIntegration &) except + nogil  # compiler

        # TODO iterator
        # TEMPLATE # void transform(InputPeakIterator begin_input, InputPeakIterator end_input, float resolution, unsigned int zeros) except + nogil 

        # parent class
        #void init(double scale, double spacing) except + nogil 

