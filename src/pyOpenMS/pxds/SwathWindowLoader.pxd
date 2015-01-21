from Types cimport *
from SwathMap cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>" namespace "OpenMS":

    cdef cppclass SwathWindowLoader:

        SwathWindowLoader() nogil except +
        SwathWindowLoader(SwathWindowLoader) nogil except + 

        void annotateSwathMapsFromFile(libcpp_string filename,
                                       libcpp_vector[ SwathMap ] & swath_maps, bool doSort) nogil except +

        void readSwathWindows(libcpp_string filename, 
                              libcpp_vector[double] & swath_prec_lower_,
                              libcpp_vector[double] & swath_prec_upper_ ) nogil except +

