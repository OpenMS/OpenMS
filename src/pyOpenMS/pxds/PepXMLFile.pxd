from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/PepXMLFile.h>" namespace "OpenMS":

    cdef cppclass PepXMLFile:

        PepXMLFile() nogil except +

        # Since PepXML may not store the complete information, it may be
        # neessary to also pass an experiment file name and an experiment MS
        # run from which retention times can be extracted.
        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String experiment_name,
                  MSExperiment[Peak1D, ChromatogramPeak] & experiment,
                  bool use_precursor_data
                  ) nogil except +

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String experiment_name
                  ) nogil except +

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) nogil except +

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) nogil except +

