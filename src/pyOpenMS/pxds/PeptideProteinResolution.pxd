from Types cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>" namespace "OpenMS":
    
    cdef cppclass PeptideProteinResolution "OpenMS::PeptideProteinResolution":
        # wrap-doc:
            #   Resolves shared peptides based on protein scores
            #   -----
            #   Resolves connected components of the bipartite protein-peptide graph based
            #   on protein probabilities/scores and adds them as additional protein_groups
            #   to the protein identification run processed.
            #   Thereby greedily assigns shared peptides in this component uniquely to the
            #   proteins of the current @em best @em indistinguishable protein group, until
            #   every peptide is uniquely assigned. This effectively allows more peptides to
            #   be used in ProteinQuantifier at the cost of potentially additional noise in
            #   the peptides quantities.
            #   In accordance with most state-of-the-art protein inference tools, only the
            #   best hit (PSM) for a peptide ID is considered.  Probability ties are
            #   currently resolved by taking the protein with larger number of peptides
            #   -----
            #   The class could provide iterator for ConnectedComponents in the
            #   future. One could extend the graph to include all PeptideHits (not only the
            #   best). It becomes a tripartite graph with larger connected components then.
            #   Maybe extend it to work with MS1 features. Separate resolution and adding
            #   groups to output

        PeptideProteinResolution(bool statistics) nogil except +
        PeptideProteinResolution(PeptideProteinResolution &) nogil except + # compiler

        void buildGraph(ProteinIdentification & protein, 
                        libcpp_vector[ PeptideIdentification ] & peptides) nogil except +
            # wrap-doc:
                #   Initialize and store the graph (= maps), needs sorted groups for
                #   correct functionality. Therefore sorts the indist. protein groups
                #   if not skipped
                #   -----
                #   :param protein: ProteinIdentification object storing IDs and groups
                #   :param peptides: Vector of ProteinIdentifications with links to the proteins
                #   :param skip_sort: Skips sorting of groups, nothing is modified then

        void resolveGraph(ProteinIdentification & protein,
                          libcpp_vector[ PeptideIdentification ] & peptides) nogil except +
            # wrap-doc:
                #   Applies resolveConnectedComponent to every component of the graph and
                #   is able to write statistics when specified. Parameters will
                #   both be mutated in this method
                #   -----
                #   :param protein: ProteinIdentification object storing IDs and groups
                #   :param peptides: vector of ProteinIdentifications with links to the proteins

        PeptideProteinResolution_ConnectedComponent findConnectedComponent(Size & root_prot_grp) nogil except +
            # wrap-doc:
                #   Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
                #   and peptides), switching from one to the other in each step
                #   -----
                #   :param root_prot_grp: Starts the BFS at this protein group index
                #   :returns: Returns a Connected Component as set of group and peptide indices

        void resolveConnectedComponent(PeptideProteinResolution_ConnectedComponent & conn_comp,
                                       ProteinIdentification & protein,
                                       libcpp_vector[ PeptideIdentification ] &
                                       peptides) nogil except +
            # wrap-doc:
                #   Resolves connected components based on posterior probabilities and adds them
                #   as additional protein_groups to the output idXML.
                #   Thereby greedily assigns shared peptides in this component uniquely to
                #   the proteins of the current BEST INDISTINGUISHABLE protein group,
                #   ready to be used in ProteinQuantifier then.
                #   This is achieved by removing all other evidence from the input
                #   PeptideIDs and iterating until each peptide is uniquely assigned.
                #   In accordance with Fido only the best hit (PSM) for an ID is considered.
                #   Probability ties resolved by taking protein with largest number of peptides
                #   -----
                #   :param conn_comp: The component to be resolved
                #   :param protein: ProteinIdentification object storing IDs and groups
                #   :param peptides: Vector of ProteinIdentifications with links to the proteins


# COMMENT: wrap static methods
cdef extern from "<OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>" namespace "OpenMS::PeptideProteinResolution":        
        # static members
        void run(libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides) nogil except +  #wrap-attach:PeptideProteinResolution

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>" namespace "OpenMS":
    
    cdef cppclass PeptideProteinResolution_ConnectedComponent "OpenMS::ConnectedComponent":
        PeptideProteinResolution_ConnectedComponent() nogil except +
        PeptideProteinResolution_ConnectedComponent(PeptideProteinResolution_ConnectedComponent) nogil except + #wrap-ignore
        libcpp_set[ size_t ] prot_grp_indices
        libcpp_set[ size_t ] pep_indices
