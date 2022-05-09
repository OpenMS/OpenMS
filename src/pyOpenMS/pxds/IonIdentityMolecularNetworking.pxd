from Types cimport *
from String cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>" namespace "OpenMS":
    cdef cppclass IonIdentityMolecularNetworking "OpenMS::IonIdentityMolecularNetworking":
        # wrap-doc:
        #  Includes the necessary functions to generate filed required for GNPS ion identity molecular networking (IIMN).

        IonIdentityMolecularNetworking() nogil except +

# wrap static methods:
cdef extern from "<OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>" namespace "OpenMS::IonIdentityMolecularNetworking":

        void annotateConsensusMap(ConsensusMap& consensus_map) nogil except + #wrap-attach:IonIdentityMolecularNetworking
        # wrap-doc:
        #  Annotate ConsensusMap for ion identity molecular networking (IIMN) workflow by GNPS.
        #  
        #  Adds meta values Constants::UserParams::IIMN_ROW_ID (unique index for each feature), Constants::UserParams::IIMN_ADDUCT_PARTNERS (related features row IDs)
        #  and Constants::UserParams::IIMN_ANNOTATION_NETWORK_NUMBER (all related features with different adduct states) get the same network number).
        #  This method requires the features annotated with the Constants::UserParams::IIMN_LINKED_GROUPS meta value.
        #  If at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN.
        #  
        #  Parameters
        #  ----------
        #  consensus_map : ConsensusMap
        #    Input ConsensusMap without IIMN annotations.

        void writeFeatureQuantificationTable(const ConsensusMap& consensus_map, const String& output_file) nogil except + #wrap-attach:IonIdentityMolecularNetworking
        # wrap-doc:
        #  Write feature quantification table (txt file) from a ConsensusMap. Required for GNPS FBMN.
        #  
        #  The table contains map information on the featureXML files from which the ConsensusMap was generated as well as
        #  a row for every consensus feature with information on rt, mz, intensity, width and quality. The same information is
        #  added for each original feature in the consensus feature.
        #  
        #  Parameters
        #  ----------
        #  consensus_map : ConsensusMap
        #    Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        #  
        #  output_file : str
        #    Output file path for the feature quantification table.

        void writeSupplementaryPairTable(const ConsensusMap& consensus_map, const String& output_file) nogil except + #wrap-attach:IonIdentityMolecularNetworking
        # wrap-doc:
        #  Write supplementary pair table (csv file) from a ConsensusMap with edge annotations for connected features. Required for GNPS IIMN.
        #  
        #  The table contains the columns "ID 1" (row ID of first feature), "ID 2" (row ID of second feature), "EdgeType" (MS1/2 annotation),
        #  "Score" (the number of direct partners from both connected features) and "Annotation" (adducts and delta m/z between two connected features).
        #  
        #  Parameters
        #  ----------
        #  consensus_map : ConsensusMap
        #    Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        #  
        #  output_file : str
        #    Output file path for the supplementary pair table.