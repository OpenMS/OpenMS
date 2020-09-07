from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>" namespace "OpenMS":
    
    cdef cppclass PeakPickerMRM(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PeakPickerMRM() nogil except +
        PeakPickerMRM(PeakPickerMRM) nogil except + #wrap-ignore

        void pickChromatogram(MSChromatogram & chromatogram, MSChromatogram & picked_chrom) nogil except +

