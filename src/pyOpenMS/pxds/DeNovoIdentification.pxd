from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoIdentification.h>" namespace "OpenMS":
    
    cdef cppclass DeNovoIdentification "OpenMS::DeNovoIdentification":
        # wrap-ignore
        # ABSTRACT CLASS
        # no-pxd-import
        DeNovoIdentification() except + nogil 
        DeNovoIdentification(DeNovoIdentification &) except + nogil 
        void getIdentifications(libcpp_vector[ PeptideIdentification ] & ids, MSExperiment & exp) except + nogil 
        void getIdentification(PeptideIdentification & id_, MSSpectrum & spectrum) except + nogil 
