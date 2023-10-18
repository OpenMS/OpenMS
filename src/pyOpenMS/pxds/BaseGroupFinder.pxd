from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>" namespace "OpenMS":
    
    cdef cppclass BaseGroupFinder(DefaultParamHandler,ProgressLogger) :
        # wrap-ignore
        # no-pxd-import
        # ABSTRACT class
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        BaseGroupFinder() except + nogil 
        BaseGroupFinder(BaseGroupFinder &) except + nogil  # compiler
        # void run(libcpp_vector[ ConsensusMap ] & input_, ConsensusMap & result) except + nogil 
        void registerChildren() except + nogil  # wrap-doc:Register all derived classes here

