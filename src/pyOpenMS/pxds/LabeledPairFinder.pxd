from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from BaseGroupFinder cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>" namespace "OpenMS":
    
    cdef cppclass LabeledPairFinder(BaseGroupFinder) :
        # wrap-inherits:
        #  BaseGroupFinder
        LabeledPairFinder() nogil except +
        LabeledPairFinder(LabeledPairFinder) nogil except + #wrap-ignore
        void run(libcpp_vector[ ConsensusMap ] & input_maps, ConsensusMap & result_map) nogil except +
        # POINTER # BaseGroupFinder * create() nogil except +
        String getProductName() nogil except +

