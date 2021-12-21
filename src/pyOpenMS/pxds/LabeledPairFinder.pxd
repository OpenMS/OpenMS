from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from BaseGroupFinder cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>" namespace "OpenMS":
    
    cdef cppclass LabeledPairFinder(BaseGroupFinder) :
        # wrap-doc:
                #   The LabeledPairFinder allows the matching of labeled features (features with a fixed distance)
                #   -----
                #   Finds feature pairs that have a defined distance in RT and m/z in the same map

        # wrap-inherits:
        #  BaseGroupFinder
        LabeledPairFinder() nogil except +
        # private
        LabeledPairFinder(LabeledPairFinder &) nogil except + # wrap-ignore
        void run(libcpp_vector[ ConsensusMap ] & input_maps, ConsensusMap & result_map) nogil except + # wrap-doc:Runs the LabeledPairFinder algorithm
        # POINTER # BaseGroupFinder * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Returns the name of this module
