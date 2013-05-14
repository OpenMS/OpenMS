from Types cimport *
from libcpp.set cimport set as libcpp_set
from MSExperiment cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>" namespace "OpenMS":
    
    cdef cppclass PeakWidthEstimator "OpenMS::PeakWidthEstimator":
        PeakWidthEstimator(PeakWidthEstimator) nogil except + #wrap-ignore
        # NAMESPACE # void estimateSpectrumFWHM(MSSpectrum[ Peak1D ] & input_, libcpp_set[ boost::tuple[ DoubleReal, DoubleReal, DoubleReal ] ] & fwhms) nogil except +
        Result estimateFWHM(MSExperiment[ Peak1D, ChromatogramPeak ] & input_) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>" namespace "OpenMS::PeakWidthEstimator":
    
    cdef cppclass Result "OpenMS::PeakWidthEstimator::Result":
        Result() nogil except +
        Result(Result) nogil except + #wrap-ignore
        DoubleReal c0
        DoubleReal c1
        Result(DoubleReal c0, DoubleReal c1) nogil except +
        # DoubleReal operator()(DoubleReal mz) nogil except +

