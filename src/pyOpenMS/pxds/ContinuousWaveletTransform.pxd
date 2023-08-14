from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Peak1D cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>" namespace "OpenMS":
    
    cdef cppclass ContinuousWaveletTransform "OpenMS::ContinuousWaveletTransform":
        ContinuousWaveletTransform() except + nogil 

        ContinuousWaveletTransform(ContinuousWaveletTransform &) except + nogil  # compiler
        libcpp_vector[ Peak1D ]  getSignal() except + nogil  # wrap-doc:Returns the wavelet transform of the signal
        void setSignal(libcpp_vector[ Peak1D ] & signal) except + nogil  # wrap-doc:Sets the wavelet transform of the signal
        libcpp_vector[ double ]  getWavelet() except + nogil  # wrap-doc:Returns the wavelet
        void setWavelet(libcpp_vector[ double ] & wavelet) except + nogil  # wrap-doc:Sets the signal
        double  getScale() except + nogil  # wrap-doc:Returns the scale of the wavelet
        void setScale(double scale) except + nogil  # wrap-doc:Sets the spacing of raw data
        double  getSpacing() except + nogil  # wrap-doc:Returns the spacing of raw data
        void setSpacing(double spacing) except + nogil  # wrap-doc:Sets the spacing of raw data
        SignedSize getLeftPaddingIndex() except + nogil  # wrap-doc:Returns the position where the signal starts (in the interval [0,end_left_padding_) are the padded zeros)
        void setLeftPaddingIndex(SignedSize end_left_padding) except + nogil  # wrap-doc:Sets the position where the signal starts
        SignedSize getRightPaddingIndex() except + nogil  # wrap-doc:Returns the position where the signal ends (in the interval (begin_right_padding_,end] are the padded zeros)
        void setRightPaddingIndex(SignedSize begin_right_padding) except + nogil  # wrap-doc:Sets the position where the signal starts
        SignedSize getSignalLength() except + nogil  # wrap-doc:Returns the signal length [end_left_padding,begin_right_padding]
        void setSignalLength(SignedSize signal_length) except + nogil  # wrap-doc:Sets the signal length [end_left_padding,begin_right_padding]
        int getSize() except + nogil  # wrap-doc:Returns the signal length including padded zeros [0,end]
        void init(double scale, double spacing) except + nogil  # wrap-doc:Perform possibly necessary preprocessing steps, like tabulating the Wavelet

        # double <](unsigned int i) except + nogil 
        # double <](unsigned int i) except + nogil 
