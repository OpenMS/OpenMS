from Types cimport *
from LinearResampler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>" namespace "OpenMS":
    
    cdef cppclass LinearResamplerAlign(LinearResampler) :
        # wrap-inherits:
        #  LinearResampler
        LinearResamplerAlign(LinearResamplerAlign) nogil except + #wrap-ignore
        # TEMPLATE # void raster(SpecT[ PeakType ] & spectrum) nogil except +
        # TEMPLATE # void raster_align(SpecT[ PeakType ] & spectrum, double start_pos, double end_pos) nogil except +
        # TEMPLATE # void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resample_end) nogil except +
        # TEMPLATE # void raster_interpolate(PeakTypeIterator raw_it, PeakTypeIterator raw_end, PeakTypeIterator it, PeakTypeIterator resampled_end) nogil except +

