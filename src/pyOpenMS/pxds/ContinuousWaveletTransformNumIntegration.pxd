from Types cimport *
# from MathFunctions cimport *
from ContinuousWaveletTransform cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>" namespace "OpenMS":
    
    cdef cppclass ContinuousWaveletTransformNumIntegration(ContinuousWaveletTransform) :
        # wrap-inherits:
        #  ContinuousWaveletTransform
        ContinuousWaveletTransformNumIntegration() nogil except +
        ContinuousWaveletTransformNumIntegration(ContinuousWaveletTransformNumIntegration) nogil except + #wrap-ignore
        # TODO iterator
        # TEMPLATE # void transform(InputPeakIterator begin_input, InputPeakIterator end_input, float resolution, unsigned int zeros) nogil except +

        # parent class
        #void init(double scale, double spacing) nogil except +

