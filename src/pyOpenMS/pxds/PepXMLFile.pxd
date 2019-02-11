from Types cimport *

from SpectrumMetaDataLookup cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FORMAT/PepXMLFile.h>" namespace "OpenMS":

    cdef cppclass PepXMLFile:

        PepXMLFile() nogil except +

        # Since PepXML may not store the complete information, it may be
        # necessary to also pass a lookup structure from which retention times
        # can be extracted.

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) nogil except +

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String experiment_name
                  ) nogil except +

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String experiment_name,
                  SpectrumMetaDataLookup lookup
                  ) nogil except +

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) nogil except +

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String mz_file,
                  String mz_name,
                  bool peptideprophet_analyzed,
                  double rt_tolerance
                  ) nogil except +

        void keepNativeSpectrumName(bool keep) nogil except +

