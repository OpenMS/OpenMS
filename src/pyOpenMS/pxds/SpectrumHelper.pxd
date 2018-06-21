from MSChromatogram cimport *
from MSSpectrum cimport *

cdef extern from "OpenMS/KERNEL/SpectrumHelper.h":

    ## Namespace only
    cdef cppclass SpectrumHelper:
        # wrap-manual-memory
        SpectrumHelper() # wrap-ignore

    void slicePeakContainer(MSChromatogram& p, double pos_start, double pos_end) # wrap-attach:SpectrumHelper
    void slicePeakContainer(MSSpectrum& p, double pos_start, double pos_end) # wrap-attach:SpectrumHelper
    void reZeroIntensities(MSChromatogram& p) # wrap-attach:SpectrumHelper
    void reZeroIntensities(MSSpectrum& p) # wrap-attach:SpectrumHelper
