from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>" namespace "OpenMS":
    cdef cppclass IonIdentityMolecularNetworking "OpenMS::IonIdentityMolecularNetworking":
        # wrap-doc:
        #  Includes the necessary functions to generate filed required for GNPS ion identity molecular networking (IIMN).

        IonIdentityMolecularNetworking() except + nogil 

# wrap static methods:
cdef extern from "<OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>" namespace "OpenMS::IonIdentityMolecularNetworking":

        void annotateConsensusMap(ConsensusMap& consensus_map) except + nogil  #wrap-attach:IonIdentityMolecularNetworking
        # wrap-doc:
        #  Annotate ConsensusMap for ion identity molecular networking (IIMN) workflow by GNPS.
        #  
        #  Adds meta values Constants::UserParams::IIMN_ROW_ID (unique index for each feature), Constants::UserParams::IIMN_ADDUCT_PARTNERS (related features row IDs)
        #  and Constants::UserParams::IIMN_ANNOTATION_NETWORK_NUMBER (all related features with different adduct states) get the same network number).
        #  This method requires the features annotated with the Constants::UserParams::IIMN_LINKED_GROUPS meta value.
        #  If at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN.
        #  
        #  
        #  :param consensus_map: Input ConsensusMap without IIMN annotations.

        void writeSupplementaryPairTable(const ConsensusMap& consensus_map, const String& output_file) except + nogil  #wrap-attach:IonIdentityMolecularNetworking
        # wrap-doc:
        #  Write supplementary pair table (csv file) from a ConsensusMap with edge annotations for connected features. Required for GNPS IIMN.
        #  
        #  The table contains the columns "ID 1" (row ID of first feature), "ID 2" (row ID of second feature), "EdgeType" (MS1/2 annotation),
        #  "Score" (the number of direct partners from both connected features) and "Annotation" (adducts and delta m/z between two connected features).
        #  
        #  
        #  :param consensus_map: Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        #  :param output_file: Output file path for the supplementary pair table.