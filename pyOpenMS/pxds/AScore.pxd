from ConsensusMap cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from PeptideHit cimport *

from MSExperiment cimport *
from MSSpectrum cimport *
from RichPeak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AScore.h>" namespace "OpenMS":

    cdef cppclass AScore:

        AScore() nogil except +
        AScore(AScore) nogil except +   # wrap-ignore

        # TODO cannot input RichPeakSpectrum...
        # void compute(PeptideHit & hit, MSSpectrum[RichPeak1D] & real_spectrum,
        #               DoubleReal fmt, int n_sites)


