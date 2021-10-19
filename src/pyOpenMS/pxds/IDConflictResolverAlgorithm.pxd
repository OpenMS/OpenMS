from Feature cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>" namespace "OpenMS":

    cdef cppclass IDConflictResolverAlgorithm(DefaultParamHandler):
        IDConflictResolverAlgorithm() nogil except + # wrap-doc:Resolves ambiguous annotations of features with peptide identifications

        IDConflictResolverAlgorithm(IDConflictResolverAlgorithm &) nogil except +

        void resolve(FeatureMap& features) nogil except +
            # wrap-doc:
            #   Resolves ambiguous annotations of features with peptide identifications
            #   -----
            #   The the filtered identifications are added to the vector of unassigned peptides
            #   and also reduced to a single best hit
            #   -----
            #   :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)

        void resolve(ConsensusMap& features) nogil except +
            # wrap-doc:
            #   Resolves ambiguous annotations of consensus features with peptide identifications
            #   -----
            #   The the filtered identifications are added to the vector of unassigned peptides
            #   and also reduced to a single best hit
            #   -----
            #   :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)

        void resolveBetweenFeatures(FeatureMap& features) nogil except +
            # wrap-doc:
            #   In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear
            #   -----
            #   This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
            #   combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
            #   appears only once, i.e. no multiplicities

        void resolveBetweenFeatures(ConsensusMap& features) nogil except +
            # wrap-doc:
            #   In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear
            #   -----
            #   This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
            #   combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
            #   appears only once, i.e. no multiplicities
