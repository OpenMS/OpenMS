# from smart_ptr cimport shared_ptr

from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *

# cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>" namespace "OpenSwath":
# 
#   # This is an abstract base class in C++ with only pure virtual functions -> we can thus not create an instance of it
#   # cdef cppclass ISpectrumAccess:
#   #     pass
#   #       ISpectrumAccess() nogil except +
#   #       ISpectrumAccess(ISpectrumAccess) nogil except +
#   #       shared_ptr[Spectrum] getSpectrumById(int id) nogil except + # wrap-ignore
#   #     
#   # ctypedef shared_ptr[ISpectrumAccess] SpectrumAccessPtr
# 
