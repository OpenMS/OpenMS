from Types cimport *
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from ISpectrumAccess cimport *

# ctypedef shared_ptr[ISpectrumAccess] SpectrumAccessPtr

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAccessTransforming(ISpectrumAccess) :
        # wrap-inherits:
        #  ISpectrumAccess
        # wrap-ignore
        # ABSTRACT class

        SpectrumAccessTransforming(SpectrumAccessTransforming) nogil except + # wrap-ignore

