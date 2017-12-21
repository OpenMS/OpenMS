from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from XMLFile cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from XQuestResultXMLFile cimport *
from CrossLinkSpectrumMatch cimport *
from PreprocessedPairSpectra cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/XQuestResultXMLFile.h>" namespace "OpenMS":
    
    cdef cppclass XQuestResultXMLFile(XMLFile) :
        # wrap-inherits:
        #  XMLFile
        XQuestResultXMLFile() nogil except +
        XQuestResultXMLFile(XQuestResultXMLFile) nogil except + #wrap-ignore

        void load(String & filename, 
                libcpp_vector[ libcpp_vector[ PeptideIdentification ] ] & csms,
                libcpp_vector[ ProteinIdentification ] & prot_ids,
                Size min_n_hits_per_spectrum,
                bool load_to_peptideHit) nogil except +

        int getNumberOfHits() nogil except +
        double getMinScore() nogil except +
        double getMaxScore() nogil except +

        void writeXQuestXML(String out_file, String base_name, 
                            libcpp_vector[ PeptideIdentification ] & peptide_ids,
                            libcpp_vector[ libcpp_vector[ OPXL_CrossLinkSpectrumMatch ] ] & all_top_csms, 
                            MSExperiment & spectra, String precursor_mass_tolerance_unit, 
                            String fragment_mass_tolerance_unit, 
                            double precursor_mass_tolerance, 
                            double fragment_mass_tolerance,
                            double fragment_mass_tolerance_xlinks, 
                            String cross_link_name, 
                            double cross_link_mass_light, 
                            DoubleList cross_link_mass_mono_link, 
                            String in_fasta, 
                            String in_decoy_fasta, 
                            StringList cross_link_residue1, 
                            StringList cross_link_residue2, 
                            double cross_link_mass_iso_shift, 
                            String enzyme_name,
                            Size missed_cleavages) nogil except +

        void writeXQuestXMLSpec(String out_file, String base_name,
                                OPXL_PreprocessedPairSpectra & preprocessed_pair_spectra, 
                                libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & spectrum_pairs, 
                                libcpp_vector[ libcpp_vector[ OPXL_CrossLinkSpectrumMatch ] ] & all_top_csms, 
                                MSExperiment & spectra) nogil except +

        void writeXQuestXMLSpec(String out_file, String base_name, 
                                libcpp_vector[ libcpp_vector[ OPXL_CrossLinkSpectrumMatch] ] & all_top_csms,
                                MSExperiment & spectra) nogil except +

