from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from PILISModel cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PILISIdentification.h>" namespace "OpenMS":
    
    cdef cppclass PILISIdentification(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PILISIdentification() nogil except +
        PILISIdentification(PILISIdentification) nogil except +
        # POINTER # void setModel(PILISModel * hmm_model) nogil except +
        # NESTED # void getIdentifications(libcpp_vector[ libcpp_map[ String, UInt ] ] & candidates, libcpp_vector[ PeptideIdentification ] & ids, MSExperiment[RichPeak1D, ChromatogramPeak] & exp) nogil except +
        # TODO # void getIdentification(libcpp_map[ String, unsigned int ] & candidates, PeptideIdentification & id_, MSSpectrum[RichPeak1D] & spectrum) nogil except +

