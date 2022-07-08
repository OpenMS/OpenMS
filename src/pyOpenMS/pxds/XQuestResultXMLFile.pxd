from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from XMLFile cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from CrossLinkSpectrumMatch cimport *
from PreprocessedPairSpectra cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/XQuestResultXMLFile.h>" namespace "OpenMS":

    cdef cppclass XQuestResultXMLFile(XMLFile) :
        # wrap-inherits:
        #  XMLFile

        XQuestResultXMLFile() nogil except +
        XQuestResultXMLFile(XQuestResultXMLFile &) nogil except +

        void load(const String & filename,
                libcpp_vector[ PeptideIdentification ] & pep_ids,
                libcpp_vector[ ProteinIdentification ] & prot_ids) nogil except +
        # wrap-doc:
                #   Load the content of the xquest.xml file into the provided data structures
                #   -----
                #   :param filename: Filename of the file which is to be loaded
                #   :param pep_ids: Where the spectra with identifications of the input file will be loaded to
                #   :param prot_ids: Where the protein identification of the input file will be loaded to

        void store(const String & filename,
                libcpp_vector[ ProteinIdentification ] & poid,
                libcpp_vector[ PeptideIdentification ] & peid) nogil except + # wrap-doc:Stores the identifications in a xQuest XML file

        int getNumberOfHits() nogil except + # wrap-doc:Returns the total number of hits in the file
        double getMinScore() nogil except + # wrap-doc:Returns minimum score among the hits in the file
        double getMaxScore() nogil except + # wrap-doc:Returns maximum score among the hits in the file

        void writeXQuestXMLSpec(const String& out_file, const String& base_name,
                                OPXL_PreprocessedPairSpectra preprocessed_pair_spectra,
                                libcpp_vector[ libcpp_pair[ size_t, size_t ] ] spectrum_pairs,
                                libcpp_vector[ libcpp_vector[ CrossLinkSpectrumMatch ] ] all_top_csms,
                                MSExperiment spectra, const bool& test_mode) nogil except +
                # wrap-doc:
                        #   Writes spec.xml output containing matching peaks between heavy and light spectra after comparing and filtering
                        #   -----
                        #   :param out_file: Path and filename for the output file
                        #   :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
                        #   :param preprocessed_pair_spectra: The preprocessed spectra after comparing and filtering
                        #   :param spectrum_pairs: Indices of spectrum pairs in the input map
                        #   :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
                        #   :param spectra: The spectra, that were searched as a PeakMap. The indices in spectrum_pairs correspond to spectra in this map

        void writeXQuestXMLSpec(const String& out_file, const String& base_name,
                                libcpp_vector[ libcpp_vector[ CrossLinkSpectrumMatch] ] all_top_csms,
                                MSExperiment spectra, const bool& test_mode) nogil except +
                # wrap-doc:
                        #   Writes spec.xml output containing spectra for visualization. This version of the function is meant to be used for label-free linkers
                        #   -----
                        #   :param out_file: Path and filename for the output file
                        #   :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
                        #   :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
                        #   :param spectra: The spectra, that were searched as a PeakMap
