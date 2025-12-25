from MSChromatogram cimport *
from MSSpectrum cimport *

cdef extern from "OpenMS/KERNEL/SpectrumHelper.h":

    ## Namespace only
    cdef cppclass SpectrumHelper:
        # wrap-manual-memory
        SpectrumHelper() # wrap-ignore

        @staticmethod
        void removePeaks(MSChromatogram& p, double pos_start, double pos_end)

        @staticmethod
        void removePeaks(MSSpectrum& p, double pos_start, double pos_end)

        @staticmethod
        void subtractMinimumIntensity(MSChromatogram& p)

        @staticmethod
        void subtractMinimumIntensity(MSSpectrum& p)
