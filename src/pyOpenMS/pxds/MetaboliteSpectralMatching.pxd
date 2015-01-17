from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from MSExperiment cimport *
from MzTab cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>" namespace "OpenMS":

    cdef cppclass MetaboliteSpectralMatching(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #    ProgressLogger
        #    DefaultParamHandler

        MetaboliteSpectralMatching() nogil except +
        MetaboliteSpectralMatching(MetaboliteSpectralMatching) nogil except + 

        double computeHyperScore(MSSpectrum[Peak1D], MSSpectrum[Peak1D], double, double) nogil except +

        void run(MSExperiment[Peak1D,ChromatogramPeak] & exp, MzTab& mz_tab) nogil except +

