from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from BaseGroupFinder cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>" namespace "OpenMS":

    cdef cppclass StablePairFinder(BaseGroupFinder) :
        # wrap-inherits:
        #  BaseGroupFinder
        StablePairFinder() except + nogil 
        # copy constructor of 'StablePairFinder' is implicitly deleted because base class 'OpenMS::BaseGroupFinder' has an inaccessible copy constructor public BaseGroupFinder
        StablePairFinder(StablePairFinder &) except + nogil  # wrap-ignore
        void run(libcpp_vector[ ConsensusMap ] & input_maps, ConsensusMap & result_map) except + nogil 
        # POINTER # BaseGroupFinder * create() except + nogil 
        
