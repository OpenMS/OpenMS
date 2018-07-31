from Peak1D cimport *
from ChromatogramPeak cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramTools.h>" namespace "OpenMS":

    cdef cppclass ChromatogramTools:
        ChromatogramTools() nogil except +

        void convertChromatogramsToSpectra(
                MSExperiment & epx
                ) nogil except +

        void convertSpectraToChromatograms(
                MSExperiment & epx,
                bool remove_spectra,
                bool force_conversion
                ) nogil except +
