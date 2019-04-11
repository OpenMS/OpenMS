from Types cimport *
from SwathMap cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>" namespace "OpenMS":

    cdef cppclass SwathWindowLoader:

        SwathWindowLoader() nogil except +
        SwathWindowLoader(SwathWindowLoader) nogil except +

        void annotateSwathMapsFromFile(String filename,
                                       libcpp_vector[ SwathMap ]& swath_maps, bool do_sort, bool force) nogil except +

        void readSwathWindows(String filename,
                              libcpp_vector[double]& swath_prec_lower,
                              libcpp_vector[double]& swath_prec_upper) nogil except +

