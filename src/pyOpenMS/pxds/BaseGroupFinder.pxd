from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>" namespace "OpenMS":
    
    cdef cppclass BaseGroupFinder(DefaultParamHandler) :
        # wrap-ignore
        # no-pxd-import
        # ABSTRACT class
        # wrap-inherits:
        #  DefaultParamHandler
        BaseGroupFinder() except + nogil 
        BaseGroupFinder(BaseGroupFinder &) except + nogil  # compiler
 
