from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/IdXMLFile.h>" namespace "OpenMS":

    cdef cppclass IdXMLFile:

        IdXMLFile() except + nogil  # wrap-doc:Used to load and store idXML files

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  PeptideIdentificationList & peptide_ids
                  ) except + nogil  # wrap-ignore

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  PeptideIdentificationList & peptide_ids,
                  String & document_id) except + nogil  # wrap-ignore

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  PeptideIdentificationList & peptide_ids,
                  String document_id) except + nogil  # wrap-ignore

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  PeptideIdentificationList & peptide_ids
                  ) except + nogil  # wrap-ignore
