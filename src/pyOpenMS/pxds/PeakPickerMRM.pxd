from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>" namespace "OpenMS":
    
    cdef cppclass PeakPickerMRM(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #  The PeakPickerMRM finds peaks a single chromatogram
        #  
        #  It uses the PeakPickerHiRes internally to find interesting seed candidates.
        #  These candidates are then expanded and a right/left border of the peak is
        #  searched
        #  Additionally, overlapping peaks can be removed
        
        PeakPickerMRM() except + nogil 
        PeakPickerMRM(PeakPickerMRM &) except + nogil 

        void pickChromatogram(MSChromatogram & chromatogram, MSChromatogram & picked_chrom) except + nogil 
            # wrap-doc:
                #  Finds peaks in a single chromatogram and annotates left/right borders
                #  
                #  It uses a modified algorithm of the PeakPickerHiRes
                #  
                #  This function will return a picked chromatogram
