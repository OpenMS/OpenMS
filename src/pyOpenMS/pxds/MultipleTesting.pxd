from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/MATH/STATISTICS/MultipleTesting.h>" namespace "OpenMS::Math":
    
    cdef cppclass MultipleTesting:
        # wrap-ignore
        # Dummy class to anchor the MultipleTesting.pyx addon module
        # The actual functionality is provided via Python functions in the addon
        MultipleTesting() except + nogil  # wrap-ignore
