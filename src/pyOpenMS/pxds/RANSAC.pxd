from Types cimport *

cdef extern from "<OpenMS/MATH/MISC/RANSAC.h>" namespace "OpenMS::Math":

    cdef cppclass RANSAC:
        pass

cdef extern from "<OpenMS/MATH/MISC/RANSAC.h>" namespace "OpenMS::Math::RANSAC":


    libcpp_vector[libcpp_pair[double,double]] ransac(
            libcpp_vector[libcpp_pair[double,double]] pairs,
            size_t n, size_t k, double t, size_t d, bool test
            ) nogil except + # wrap-attach:RANSAC

    double llsm_rsq(
            libcpp_vector[libcpp_pair[double,double]] pairs
            ) nogil except + # wrap-attach:RANSAC


