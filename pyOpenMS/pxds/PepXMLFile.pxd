from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FORMAT/PepXMLFile.h>" namespace "OpenMS":

    cdef cppclass PepXMLFile:

        PepXMLFile() nogil except +

        void load(String & filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String & experiment_name,
                  MSExperiment[Peak1D, ChromatogramPeak] & experiment,
                  bool use_precursor_data
                  )

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String & experiment_name
                  )

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  )

