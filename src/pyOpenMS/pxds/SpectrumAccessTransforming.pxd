from Types cimport *
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from ISpectrumAccess cimport *

# ctypedef shared_ptr[ISpectrumAccess] SpectrumAccessPtr

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAccessTransforming(ISpectrumAccess) :
        # wrap-inherits:
        #  ISpectrumAccess

        #TODO: This was included in wrap-ignore with autowrap 22. 22 just didn't ignore it. 23 does, so we remove the annotiation since derived classes refer to it.


        SpectrumAccessTransforming(SpectrumAccessTransforming) except + nogil  # wrap-ignore

