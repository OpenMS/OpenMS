from Types cimport *
# from RANSACModelLinear cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSAC.h>" namespace "OpenMS::Math":

    cdef cppclass RANSACParam:

        RANSACParam() nogil except +
        RANSACParam(RANSACParam) nogil except + #wrap-ignore
        RANSACParam(size_t p_n, size_t p_k, double p_t, size_t p_d, bool p_relative_d) nogil except +

        libcpp_string toString() nogil except +

        size_t n #; //< data points; The minimum number of data points required to fit the model
        size_t k # ; //< iterations; The maximum number of iterations allowed in the algorithm 
        double t # ; //< Threshold value; for determining when a data point fits a model. Corresponds to the maximal squared deviation in units of the _second_ dimension (dim2).
        size_t d # ; //< The number of close data values (according to 't') required to assert that a model fits well to data
        bool relative_d #; //< Should 'd' be interpreted as percentages (0-100) of data input size.
        # int (*rng)(int); //< Optional RNG function (useful for testing with fixed seeds)

