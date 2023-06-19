from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from Peak1D cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>" namespace "OpenMS":
    
    cdef cppclass MRMFragmentSelection(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MRMFragmentSelection() nogil except +
        MRMFragmentSelection(MRMFragmentSelection&) nogil except +
        void selectFragments(libcpp_vector[ Peak1D ] & selected_peaks, MSSpectrum & spec) nogil except + # wrap-doc:Selects accordingly to the parameters the best peaks of spec and writes them into `selected_peaks`


