from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>" namespace "OpenMS":
    
    cdef cppclass PeakPickerChromatogram(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #  The PeakPickerChromatogram finds peaks a single chromatogram
        #  
        #  It uses the PeakPickerHiRes internally to find interesting seed candidates.
        #  These candidates are then expanded and a right/left border of the peak is
        #  searched
        #  Additionally, overlapping peaks can be removed
        
        PeakPickerChromatogram() except + nogil 
        PeakPickerChromatogram(PeakPickerChromatogram &) except + nogil 

        void pickChromatogram(MSChromatogram & chromatogram, MSChromatogram & picked_chrom) except + nogil 
            # wrap-doc:
                #  Finds peaks in a single chromatogram and annotates left/right borders
                #  
                #  It uses a modified algorithm of the PeakPickerHiRes
                #  
                #  This function will return a picked chromatogram
