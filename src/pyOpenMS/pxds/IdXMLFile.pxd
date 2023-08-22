from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/IdXMLFile.h>" namespace "OpenMS":

    cdef cppclass IdXMLFile:

        IdXMLFile() except + nogil  # wrap-doc:Used to load and store idXML files

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  ) except + nogil  # wrap-doc:Loads the identifications of an idXML file without identifier

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String document_id) except + nogil  # wrap-doc:Stores the data in an idXML file

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids
                  ) except + nogil 
