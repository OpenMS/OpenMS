from Peak1D cimport *
from ChromatogramPeak cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramTools.h>" namespace "OpenMS":

    cdef cppclass ChromatogramTools:
        ChromatogramTools() nogil except +
        ChromatogramTools(ChromatogramTools &) nogil except +

        void convertChromatogramsToSpectra(
                MSExperiment & epx
                ) nogil except + # wrap-doc:Converts the chromatogram to a list of spectra with instrument settings

        void convertSpectraToChromatograms(
                MSExperiment & epx,
                bool remove_spectra,
                bool force_conversion
                ) nogil except + # wrap-doc:Converts e.g. SRM spectra to chromatograms
