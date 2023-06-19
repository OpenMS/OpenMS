from Types cimport *
from MultiplexDeltaMasses cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexIsotopicPeakPattern "OpenMS::MultiplexIsotopicPeakPattern":

        MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi) nogil except +
        MultiplexIsotopicPeakPattern(MultiplexIsotopicPeakPattern &) nogil except + # compiler

        int getCharge() nogil except + # wrap-doc:Returns charge
        int getPeaksPerPeptide() nogil except + # wrap-doc:Returns peaks per peptide
        MultiplexDeltaMasses getMassShifts() nogil except + # wrap-doc:Returns mass shifts
        int getMassShiftIndex() nogil except + # wrap-doc:Returns mass shift index
        unsigned getMassShiftCount() nogil except + # wrap-doc:Returns number of mass shifts i.e. the number of peptides in the multiplet
        double getMassShiftAt(int i) nogil except + # wrap-doc:Returns mass shift at position i
        double getMZShiftAt(int i) nogil except + # wrap-doc:Returns m/z shift at position i
        unsigned getMZShiftCount() nogil except + # wrap-doc:Returns number of m/z shifts

