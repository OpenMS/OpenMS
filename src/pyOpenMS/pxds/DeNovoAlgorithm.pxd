from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from DeNovoIonScoring cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass DeNovoAlgorithm "OpenMS::DeNovoAlgorithm":
        # wrap-ignore
        DeNovoAlgorithm() nogil except +
        DeNovoAlgorithm(DeNovoAlgorithm) nogil except +
        void generateCandidates(libcpp_vector[ PeptideIdentification ] & candidates, libcpp_vector[ libcpp_vector[ IonScore_DeNovoIonScoring ] ] &ion_scores, MSExperiment[RichPeak1D, ChromatogramPeak] &exp) nogil except +
        void generateCandidates(PeptideIdentification &candidates, libcpp_vector[ IonScore_DeNovoIonScoring ] & ion_scores, MSSpectrum[RichPeak1D] & spec) nogil except +

