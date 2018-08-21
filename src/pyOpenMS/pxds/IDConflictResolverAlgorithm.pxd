from Feature cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>" namespace "OpenMS":

    cdef cppclass IDConflictResolverAlgorithm(DefaultParamHandler):
        IDConflictResolverAlgorithm() nogil except +
        
        IDConflictResolverAlgorithm(IDConflictResolverAlgorithm) nogil except + # wrap-ignore

        void resolve(FeatureMap& features) nogil except +

        void resolve(ConsensusMap& features) nogil except +

        void resolveBetweenFeatures(FeatureMap& features) nogil except +

        void resolveBetweenFeatures(ConsensusMap& features) nogil except +
