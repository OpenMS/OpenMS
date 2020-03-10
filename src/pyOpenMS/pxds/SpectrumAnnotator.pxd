from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from TheoreticalSpectrumGenerator cimport *
from SpectrumAlignment cimport *
from Peak1D cimport *
from PeptideHit cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/CHEMISTRY/SpectrumAnnotator.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAnnotator(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        SpectrumAnnotator() nogil except +
        SpectrumAnnotator(SpectrumAnnotator) nogil except +

        void annotateMatches(MSSpectrum & spec, PeptideHit & ph, TheoreticalSpectrumGenerator & tg, SpectrumAlignment & sa) nogil except +
        void addIonMatchStatistics(PeptideIdentification & pi, MSSpectrum & spec, TheoreticalSpectrumGenerator & tg, SpectrumAlignment & sa) nogil except +

