from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from RichPeak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoIonScoring.h>" namespace "OpenMS":
    
    cdef cppclass DeNovoIonScoring "OpenMS::DeNovoIonScoring":
        # wrap-ignore
        DeNovoIonScoring() nogil except +
        DeNovoIonScoring(DeNovoIonScoring) nogil except +
        void getIonScores(libcpp_vector[ IonScore_DeNovoIonScoring ] &ion_scores, MSSpectrum[RichPeak1D] &spec) nogil except +
        void getIonScores(libcpp_vector[ libcpp_vector[ IonScore_DeNovoIonScoring ] ] &ion_scores, MSExperiment[RichPeak1D, ChromatogramPeak] &exp) nogil except + #wrap-ignore


cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoIonScoring.h>" namespace "OpenMS::DeNovoIonScoring":
    
    cdef cppclass IonScore_DeNovoIonScoring "OpenMS::DeNovoIonScoring::IonScore":
        IonScore_DeNovoIonScoring() nogil except +
        IonScore_DeNovoIonScoring(IonScore_DeNovoIonScoring) nogil except +

        # score
        DoubleReal score

        # position of the ion
        DoubleReal position

        # index of peak in the spectrum, -1 if not in spectrum
        ptrdiff_t index

