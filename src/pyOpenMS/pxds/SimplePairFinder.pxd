from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from BaseGroupFinder cimport *
from ConsensusMap cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>" namespace "OpenMS":
    
    cdef cppclass SimplePairFinder(BaseGroupFinder) :
        # wrap-inherits:
        #  BaseGroupFinder
        # wrap-doc:
        # This class implements a simple point pair finding algorithm

        SimplePairFinder() nogil except +
        # copy constructor of 'SimplePairFinder' is implicitly deleted because base class 'OpenMS::BaseGroupFinder' has an inaccessible copy constructor public BaseGroupFinder
        SimplePairFinder(SimplePairFinder &) nogil except + # wrap-ignore

        void run(libcpp_vector[ ConsensusMap ] & input_maps, ConsensusMap & result_map) nogil except +
        # POINTER # BaseGroupFinder * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Returns the name of this module
