from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoIonScoring.h>" namespace "OpenMS":
    
    cdef cppclass DeNovoIonScoring "OpenMS::DeNovoIonScoring":
        # wrap-ignore
        # ABSTRACT CLASS
        # no-pxd-import
        DeNovoIonScoring() nogil except +
        DeNovoIonScoring(DeNovoIonScoring) nogil except +
        void getIonScores(libcpp_vector[ IonScore ] &ion_scores, MSSpectrum[Peak1D] &spec) nogil except +
        void getIonScores(libcpp_vector[ libcpp_vector[ IonScore ] ] &ion_scores, MSExperiment &exp) nogil except + #wrap-ignore

cdef extern from "<OpenMS/ANALYSIS/DENOVO/DeNovoIonScoring.h>" namespace "OpenMS::DeNovoIonScoring":
    
    cdef cppclass IonScore "OpenMS::DeNovoIonScoring::IonScore":
        # wrap-attach:
        #    DeNovoIonScoring
        IonScore() nogil except +
        IonScore(IonScore) nogil except +

        # score
        double score

        # position of the ion
        double position

        # index of peak in the spectrum, -1 if not in spectrum
        ptrdiff_t index
