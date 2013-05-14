from Peak1D cimport *
from ChromatogramPeak cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramTools.h>" namespace "OpenMS":

    cdef cppclass ChromatogramTools:
        ChromatogramTools()

        void convertChromatogramsToSpectra(MSExperiment[Peak1D, ChromatogramPeak] & epx) nogil except +
        void convertSpectraToChromatograms(MSExperiment[Peak1D, ChromatogramPeak] & epx, int remove_spectra) nogil except +
