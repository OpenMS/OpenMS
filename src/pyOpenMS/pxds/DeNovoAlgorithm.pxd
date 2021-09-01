from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from DeNovoIonScoring cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoAlgorithm.h>" namespace "OpenMS":

    cdef cppclass DeNovoAlgorithm "OpenMS::DeNovoAlgorithm":
        # wrap-ignore
        # no-pxd-import
        DeNovoAlgorithm() nogil except +
        DeNovoAlgorithm(DeNovoAlgorithm &) nogil except +
        void generateCandidates(libcpp_vector[ PeptideIdentification ] & candidates,
                                libcpp_vector[ libcpp_vector[ DeNovoIonScore ] ] &ion_scores,
                                MSExperiment &exp) nogil except +
        void generateCandidates(PeptideIdentification &candidates,
                                libcpp_vector[ DeNovoIonScore ] & ion_scores,
                                MSSpectrum & spec) nogil except +

