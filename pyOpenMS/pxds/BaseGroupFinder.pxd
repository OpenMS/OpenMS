from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from ConsensusMap cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>" namespace "OpenMS":
    
    cdef cppclass BaseGroupFinder(DefaultParamHandler,ProgressLogger) :
        # wrap-ignore
        # ABSTRACT class
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        BaseGroupFinder() nogil except +
        BaseGroupFinder(BaseGroupFinder) nogil except + #wrap-ignore
        # void run(libcpp_vector[ ConsensusMap ] & input_, ConsensusMap & result) nogil except +
        void registerChildren() nogil except +

