from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/IdXMLFile.h>" namespace "OpenMS":

    cdef cppclass IdXMLFile:
        # wrap-doc:
        #  File adapter for idXML files
        #  
        #  Provides methods to load and store identification data in idXML format.
        #  idXML files store protein and peptide identifications from database search engines.
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    protein_ids = []
        #    peptide_ids = []
        #    IdXMLFile().load("test.idXML", protein_ids, peptide_ids)

        IdXMLFile() except + nogil

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
