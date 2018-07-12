from MSChromatogram cimport *
from MSSpectrum cimport *

cdef extern from "OpenMS/KERNEL/SpectrumHelper.h":

    ## Namespace only
    cdef cppclass SpectrumHelper:
        # wrap-manual-memory
        SpectrumHelper() # wrap-ignore

    void removePeaks(MSChromatogram& p, double pos_start, double pos_end) # wrap-attach:SpectrumHelper
    void removePeaks(MSSpectrum& p, double pos_start, double pos_end) # wrap-attach:SpectrumHelper
    void subtractMinimumIntensity(MSChromatogram& p) # wrap-attach:SpectrumHelper
    void subtractMinimumIntensity(MSSpectrum& p) # wrap-attach:SpectrumHelper
