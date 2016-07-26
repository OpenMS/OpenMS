from Types cimport *
from MultiplexDeltaMasses cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexIsotopicPeakPattern "OpenMS::MultiplexIsotopicPeakPattern":

        MultiplexIsotopicPeakPattern(MultiplexIsotopicPeakPattern) nogil except + #wrap-ignore
        MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi) nogil except +

        int getCharge() nogil except +
        int getPeaksPerPeptide() nogil except +
        MultiplexDeltaMasses getMassShifts() nogil except +
        int getMassShiftIndex() nogil except +
        unsigned getMassShiftCount() nogil except +
        double getMassShiftAt(int i) nogil except +
        double getMZShiftAt(int i) nogil except +
        unsigned getMZShiftCount() nogil except +

