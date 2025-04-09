from Types cimport *
from String cimport *
from RANSACModelLinear cimport *
from RANSACModelQuadratic cimport *

cdef extern from "<OpenMS/ML/RANSAC/RANSAC.h>" namespace "OpenMS::Math":

    cdef cppclass RANSAC[TModelType]:
        # wrap-instances:
        #  RANSAC := RANSAC[RansacModelLinear]
        #  RANSACQuadratic := RANSAC[RansacModelQuadratic]

        RANSAC() except + nogil 
        RANSAC(UInt64 seed) except + nogil 
        RANSAC(RANSAC[TModelType] &) except + nogil  # wrap-ignore
        
        libcpp_vector[libcpp_pair[double,double]] ransac(
            libcpp_vector[libcpp_pair[double,double]] pairs,
            size_t n, size_t k, double t, size_t d, bool relative_d
            ) except + nogil 

    cdef cppclass RANSACParam:

        RANSACParam() except + nogil  # wrap-doc:A simple struct to carry all the parameters required for a RANSAC run
        RANSACParam(RANSACParam) except + nogil  # wrap-ignore
        RANSACParam(size_t p_n, size_t p_k, double p_t, size_t p_d, bool p_relative_d) except + nogil 

        String toString() except + nogil 

        size_t n #; //< data points; The minimum number of data points required to fit the model
        size_t k # ; //< iterations; The maximum number of iterations allowed in the algorithm 
        double t # ; //< Threshold value; for determining when a data point fits a model. Corresponds to the maximal squared deviation in units of the _second_ dimension (dim2).
        size_t d # ; //< The number of close data values (according to 't') required to assert that a model fits well to data
        bool relative_d #; //< Should 'd' be interpreted as percentages (0-100) of data input size.
        # int (*rng)(int); //< Optional RNG function (useful for testing with fixed seeds)

