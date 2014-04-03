from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *
from RichPeak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoIdentification.h>" namespace "OpenMS":
    
    cdef cppclass DeNovoIdentification "OpenMS::DeNovoIdentification":
        # wrap-ignore
        DeNovoIdentification() nogil except +
        DeNovoIdentification(DeNovoIdentification) nogil except +
        void getIdentifications(libcpp_vector[ PeptideIdentification ] & ids, MSExperiment[RichPeak1D, ChromatogramPeak] & exp) nogil except +
        void getIdentification(PeptideIdentification & id_, MSSpectrum[RichPeak1D] & spectrum) nogil except +

