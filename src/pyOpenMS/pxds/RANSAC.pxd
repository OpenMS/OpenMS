from Types cimport *
from String cimport *
from RANSACModelLinear cimport *
from RANSACModelQuadratic cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSAC.h>" namespace "OpenMS::Math":

    cdef cppclass RANSAC[TModelType]:
        # wrap-instances:
        #   RANSAC := RANSAC[RansacModelLinear]
        #   RANSACQuadratic := RANSAC[RansacModelQuadratic]

        RANSAC() nogil except +
        RANSAC(RANSAC[TModelType] &) nogil except + # wrap-ignore

    cdef cppclass RANSACParam:

        RANSACParam() nogil except +
        RANSACParam(RANSACParam) nogil except + # wrap-ignore
        RANSACParam(size_t p_n, size_t p_k, double p_t, size_t p_d, bool p_relative_d) nogil except +

        String toString() nogil except +

        size_t n #; //< data points; The minimum number of data points required to fit the model
        size_t k # ; //< iterations; The maximum number of iterations allowed in the algorithm 
        double t # ; //< Threshold value; for determining when a data point fits a model. Corresponds to the maximal squared deviation in units of the _second_ dimension (dim2).
        size_t d # ; //< The number of close data values (according to 't') required to assert that a model fits well to data
        bool relative_d #; //< Should 'd' be interpreted as percentages (0-100) of data input size.
        # int (*rng)(int); //< Optional RNG function (useful for testing with fixed seeds)

cdef extern from "<OpenMS/MATH/MISC/RANSAC.h>" namespace "OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>":

    libcpp_vector[libcpp_pair[double,double]] ransac(
            libcpp_vector[libcpp_pair[double,double]] pairs,
            size_t n, size_t k, double t, size_t d, bool test
            ) nogil except + # wrap-attach:RANSAC

cdef extern from "<OpenMS/MATH/MISC/RANSAC.h>" namespace "OpenMS::Math::RANSAC<OpenMS:Math::RansacModelQuadratic>":

    libcpp_vector[libcpp_pair[double,double]] ransac(
            libcpp_vector[libcpp_pair[double,double]] pairs,
            size_t n, size_t k, double t, size_t d, bool test
            ) nogil except + # wrap-attach:RANSACQuadratic

