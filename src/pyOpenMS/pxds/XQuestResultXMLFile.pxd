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
                libcpp_vector[ PeptideIdentification ] & pep_ids,
                libcpp_vector[ ProteinIdentification ] & prot_ids) nogil except +

        void store(String & filename,
                libcpp_vector[ ProteinIdentification ] & poid,
                libcpp_vector[ PeptideIdentification ] & peid) nogil except +

        int getNumberOfHits() nogil except +
        double getMinScore() nogil except +
        double getMaxScore() nogil except +

        void writeXQuestXMLSpec(String out_file, String base_name,
                                OPXL_PreprocessedPairSpectra & preprocessed_pair_spectra,
                                libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & spectrum_pairs,
                                libcpp_vector[ libcpp_vector[ OPXL_CrossLinkSpectrumMatch ] ] & all_top_csms,
                                MSExperiment & spectra) nogil except +

        void writeXQuestXMLSpec(String out_file, String base_name,
                                libcpp_vector[ libcpp_vector[ OPXL_CrossLinkSpectrumMatch] ] & all_top_csms,
                                MSExperiment & spectra) nogil except +
