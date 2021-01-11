from Types cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>" namespace "OpenMS":
    
    cdef cppclass PeptideProteinResolution "OpenMS::PeptideProteinResolution":
        PeptideProteinResolution(PeptideProteinResolution) nogil except + #wrap-ignore
        PeptideProteinResolution(bool statistics) nogil except +

        void buildGraph(ProteinIdentification & protein, 
                        libcpp_vector[ PeptideIdentification ] & peptides) nogil except +

        void resolveGraph(ProteinIdentification & protein,
                          libcpp_vector[ PeptideIdentification ] & peptides) nogil except +

        PeptideProteinResolution_ConnectedComponent findConnectedComponent(Size & root_prot_grp) nogil except +

        void resolveConnectedComponent(PeptideProteinResolution_ConnectedComponent & conn_comp,
                                       ProteinIdentification & protein,
                                       libcpp_vector[ PeptideIdentification ] &
                                       peptides) nogil except +


cdef extern from "<OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>" namespace "OpenMS":
    
    cdef cppclass PeptideProteinResolution_ConnectedComponent "OpenMS::ConnectedComponent":
        PeptideProteinResolution_ConnectedComponent() nogil except +
        PeptideProteinResolution_ConnectedComponent(PeptideProteinResolution_ConnectedComponent) nogil except + #wrap-ignore
        libcpp_set[ size_t ] prot_grp_indices
        libcpp_set[ size_t ] pep_indices
