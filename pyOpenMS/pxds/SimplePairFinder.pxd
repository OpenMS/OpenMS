from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
# from BaseGroupFinder cimport *
from ConsensusMap cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>" namespace "OpenMS":
    
    cdef cppclass SimplePairFinder(BaseGroupFinder) :
        # wrap-inherits:
        #  BaseGroupFinder
        SimplePairFinder() nogil except +
        SimplePairFinder(SimplePairFinder) nogil except + #wrap-ignore
        void run(libcpp_vector[ ConsensusMap ] & input_maps, ConsensusMap & result_map) nogil except +
        # POINTER # BaseGroupFinder * create() nogil except +
        String getProductName() nogil except +
