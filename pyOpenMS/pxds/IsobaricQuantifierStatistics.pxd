from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricQuantifierStatistics "OpenMS::IsobaricQuantifierStatistics":
        IsobaricQuantifierStatistics() nogil except +
        IsobaricQuantifierStatistics(IsobaricQuantifierStatistics) nogil except +
        Size channel_count
        Size iso_number_ms2_negative
        Size iso_number_reporter_negative
        Size iso_number_reporter_different
        DoubleReal iso_solution_different_intensity
        DoubleReal iso_total_intensity_negative
        Size number_ms2_total
        Size number_ms2_empty
        libcpp_map[ size_t, size_t ] empty_channels
        void reset() nogil except +

