from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from Peak1D cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>" namespace "OpenMS":
    
    cdef cppclass MRMFragmentSelection(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        MRMFragmentSelection() except + nogil 
        MRMFragmentSelection(MRMFragmentSelection&) except + nogil 
        void selectFragments(libcpp_vector[ Peak1D ] & selected_peaks, MSSpectrum & spec) except + nogil  # wrap-doc:Selects accordingly to the parameters the best peaks of spec and writes them into `selected_peaks`


