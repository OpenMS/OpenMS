from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>" namespace "OpenMS":
    
    cdef cppclass ContinuousWaveletTransform "OpenMS::ContinuousWaveletTransform":
        ContinuousWaveletTransform() nogil except +
        ContinuousWaveletTransform(ContinuousWaveletTransform) nogil except + #wrap-ignore
        libcpp_vector[ Peak1D ]  getSignal() nogil except +
        void setSignal(libcpp_vector[ Peak1D ] & signal) nogil except +
        libcpp_vector[ double ]  getWavelet() nogil except +
        void setWavelet(libcpp_vector[ double ] & wavelet) nogil except +
        double  getScale() nogil except +
        void setScale(DoubleReal scale) nogil except +
        double  getSpacing() nogil except +
        void setSpacing(double spacing) nogil except +
        SignedSize getLeftPaddingIndex() nogil except +
        void setLeftPaddingIndex(SignedSize end_left_padding) nogil except +
        SignedSize getRightPaddingIndex() nogil except +
        void setRightPaddingIndex(SignedSize begin_right_padding) nogil except +
        SignedSize getSignalLength() nogil except +
        void setSignalLength(SignedSize signal_length) nogil except +
        int getSize() nogil except +
        void init(double scale, double spacing) nogil except +
        # double <](unsigned int i) nogil except +
        # double <](unsigned int i) nogil except +

