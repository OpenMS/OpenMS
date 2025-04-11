from Peak1D cimport *
from ChromatogramPeak cimport *
from MSRun cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramTools.h>" namespace "OpenMS":

    cdef cppclass ChromatogramTools:
        ChromatogramTools() except + nogil 
        ChromatogramTools(ChromatogramTools &) except + nogil 

        void convertChromatogramsToSpectra(
                MSRun & epx
                ) except + nogil  # wrap-doc:Converts the chromatogram to a list of spectra with instrument settings

        void convertSpectraToChromatograms(
                MSRun & epx,
                bool remove_spectra,
                bool force_conversion
                ) except + nogil  # wrap-doc:Converts e.g. SRM spectra to chromatograms
