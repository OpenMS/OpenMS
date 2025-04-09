from Types cimport *
from LinearResampler cimport *

cdef extern from "<OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>" namespace "OpenMS":
    
    cdef cppclass LinearResamplerAlign(LinearResampler) :
        # wrap-inherits:
        #  LinearResampler
        LinearResamplerAlign(LinearResamplerAlign &) except + nogil 

        # TEMPLATE # void raster(SpecT[ PeakType ] & spectrum) except + nogil 
        # TEMPLATE # void raster_align(SpecT[ PeakType ] & spectrum, double start_pos, double end_pos) except + nogil 
        # TEMPLATE # void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resample_end) except + nogil 
        # TEMPLATE # void raster_interpolate(PeakTypeIterator raw_it, PeakTypeIterator raw_end, PeakTypeIterator it, PeakTypeIterator resampled_end) except + nogil 

