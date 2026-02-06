from Types cimport *

from String cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>" namespace "OpenMS":

    cdef cppclass IsobaricQuantifierStatistics "OpenMS::IsobaricQuantifierStatistics":
        # wrap-doc:
        #  Statistics for isobaric quantitation performance. Tracks the number of
        #  MS2 spectra processed, empty channels, and comparison metrics between
        #  NNLS (non-negative least squares) and naive matrix inversion methods

        IsobaricQuantifierStatistics() except + nogil  # wrap-doc:Default constructor
        IsobaricQuantifierStatistics(IsobaricQuantifierStatistics &) except + nogil  # compiler
        Size channel_count  # wrap-doc:Number of channels in the quantitation method (e.g., 4, 6, 8, 10, 11, 16, 18)
        Size iso_number_ms2_negative  # wrap-doc:Number of MS2 spectra where one or more channels had negative solution
        Size iso_number_reporter_negative  # wrap-doc:Number of channels where naive solution was negative
        Size iso_number_reporter_different  # wrap-doc:Number of channels greater than 0 where naive solution differed from NNLS
        double iso_solution_different_intensity  # wrap-doc:Absolute intensity difference between NNLS and naive solutions (for channels greater than 0)
        double iso_total_intensity_negative  # wrap-doc:Total intensity for spectra where naive solution produced negative values
        Size number_ms2_total  # wrap-doc:Total number of MS2 spectra processed
        Size number_ms2_empty  # wrap-doc:Number of empty MS2 spectra (no reporter ions detected)
        # libcpp_map[ size_t, size_t ] empty_channels
        void reset() except + nogil  # wrap-doc:Resets all statistics counters to zero
