from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from BaseGroupFinder cimport *
from ProgressLogger cimport *
# from HashGrid cimport *
# from GridFeature cimport *
from QTCluster cimport *
from FeatureDistance cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>" namespace "OpenMS":
    
    cdef cppclass QTClusterFinder(BaseGroupFinder) :
        # wrap-inherits:
        #  BaseGroupFinder
        QTClusterFinder() except + nogil 
        # copy constructor of 'QTClusterFinder' is implicitly deleted because base class 'OpenMS::BaseGroupFinder' has an inaccessible copy constructor public BaseGroupFinder
        QTClusterFinder(QTClusterFinder &) except + nogil  # wrap-ignore
        void run(libcpp_vector[ ConsensusMap ] & input_maps, ConsensusMap & result_map) except + nogil 
        void run(libcpp_vector[ FeatureMap ] & input_maps, ConsensusMap & result_map) except + nogil 
        # POINTER # BaseGroupFinder * create() except + nogil 
