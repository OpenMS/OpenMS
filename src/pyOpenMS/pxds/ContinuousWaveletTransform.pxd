from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>" namespace "OpenMS":
    
    cdef cppclass ContinuousWaveletTransform "OpenMS::ContinuousWaveletTransform":
        ContinuousWaveletTransform() nogil except +

        ContinuousWaveletTransform(ContinuousWaveletTransform &) nogil except + # compiler
        libcpp_vector[ Peak1D ]  getSignal() nogil except + # wrap-doc:Returns the wavelet transform of the signal
        void setSignal(libcpp_vector[ Peak1D ] & signal) nogil except + # wrap-doc:Sets the wavelet transform of the signal
        libcpp_vector[ double ]  getWavelet() nogil except + # wrap-doc:Returns the wavelet
        void setWavelet(libcpp_vector[ double ] & wavelet) nogil except + # wrap-doc:Sets the signal
        double  getScale() nogil except + # wrap-doc:Returns the scale of the wavelet
        void setScale(double scale) nogil except + # wrap-doc:Sets the spacing of raw data
        double  getSpacing() nogil except + # wrap-doc:Returns the spacing of raw data
        void setSpacing(double spacing) nogil except + # wrap-doc:Sets the spacing of raw data
        SignedSize getLeftPaddingIndex() nogil except + # wrap-doc:Returns the position where the signal starts (in the interval [0,end_left_padding_) are the padded zeros)
        void setLeftPaddingIndex(SignedSize end_left_padding) nogil except + # wrap-doc:Sets the position where the signal starts
        SignedSize getRightPaddingIndex() nogil except + # wrap-doc:Returns the position where the signal ends (in the interval (begin_right_padding_,end] are the padded zeros)
        void setRightPaddingIndex(SignedSize begin_right_padding) nogil except + # wrap-doc:Sets the position where the signal starts
        SignedSize getSignalLength() nogil except + # wrap-doc:Returns the signal length [end_left_padding,begin_right_padding]
        void setSignalLength(SignedSize signal_length) nogil except + # wrap-doc:Sets the signal length [end_left_padding,begin_right_padding]
        int getSize() nogil except + # wrap-doc:Returns the signal length including padded zeros [0,end]
        void init(double scale, double spacing) nogil except + # wrap-doc:Perform possibly necessary preprocessing steps, like tabulating the Wavelet

        # double <](unsigned int i) nogil except +
        # double <](unsigned int i) nogil except +
