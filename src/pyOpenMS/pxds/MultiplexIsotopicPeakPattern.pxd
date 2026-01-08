from Types cimport *
from MultiplexDeltaMasses cimport *

cdef extern from "<OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>" namespace "OpenMS":
    
    cdef cppclass MultiplexIsotopicPeakPattern "OpenMS::MultiplexIsotopicPeakPattern":

        MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi) except + nogil 
        MultiplexIsotopicPeakPattern(MultiplexIsotopicPeakPattern &) except + nogil  # compiler

        int getCharge() except + nogil  # wrap-doc:Returns charge
        int getPeaksPerPeptide() except + nogil  # wrap-doc:Returns peaks per peptide
        MultiplexDeltaMasses getMassShifts() except + nogil  # wrap-doc:Returns mass shifts
        int getMassShiftIndex() except + nogil  # wrap-doc:Returns mass shift index
        unsigned getMassShiftCount() except + nogil  # wrap-doc:Returns number of mass shifts i.e. the number of peptides in the multiplet
        double getMassShiftAt(int i) except + nogil  # wrap-doc:Returns mass shift at position i
        double getMZShiftAt(int i) except + nogil  # wrap-doc:Returns m/z shift at position i
        unsigned getMZShiftCount() except + nogil  # wrap-doc:Returns number of m/z shifts

