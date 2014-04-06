from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from RichPeak1D cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>" namespace "OpenMS":
    
    cdef cppclass MRMFragmentSelection(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MRMFragmentSelection() nogil except +
        MRMFragmentSelection(MRMFragmentSelection) nogil except +
        void selectFragments(libcpp_vector[ RichPeak1D ] & selected_peaks, MSSpectrum[RichPeak1D] & spec) nogil except +

