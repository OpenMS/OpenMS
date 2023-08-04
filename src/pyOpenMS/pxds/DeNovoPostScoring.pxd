from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoPostScoring.h>" namespace "OpenMS":
    
    cdef cppclass DeNovoPostScoring "OpenMS::DeNovoPostScoring":
        # wrap-ignore
        # ABSTRACT CLASS
        # no-pxd-import
        DeNovoPostScoring() except + nogil 
        DeNovoPostScoring(DeNovoPostScoring &) except + nogil 
        void apply(libcpp_vector[ PeptideIdentification ] &identifications, MSExperiment &exp) except + nogil 
        void apply(PeptideIdentification &identification, MSSpectrum &spec) except + nogil 

