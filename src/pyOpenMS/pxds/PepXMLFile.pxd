from Types cimport *

from SpectrumMetaDataLookup cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FORMAT/PepXMLFile.h>" namespace "OpenMS":

    cdef cppclass PepXMLFile:

        PepXMLFile() except + nogil 
        #  copy constructor of 'PepXMLFile' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor protected Internal::XMLHandler,
        PepXMLFile(PepXMLFile &) except + nogil  # wrap-ignore

        # Since PepXML may not store the complete information, it may be
        # necessary to also pass a lookup structure from which retention times
        # can be extracted.

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) except + nogil 

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String experiment_name
                  ) except + nogil 

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String experiment_name,
                  SpectrumMetaDataLookup lookup
                  ) except + nogil 

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) except + nogil 

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String mz_file,
                  String mz_name,
                  bool peptideprophet_analyzed,
                  double rt_tolerance
                  ) except + nogil 

        void keepNativeSpectrumName(bool keep) except + nogil 
        void setParseUnknownScores(bool parse_unknown_scores) except + nogil 

