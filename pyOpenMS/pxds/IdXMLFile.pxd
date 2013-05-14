from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/FORMAT/IdXMLFile.h>" namespace "OpenMS":

    cdef cppclass IdXMLFile:

        IdXMLFile() nogil except +

        void load(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  )

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  String document_id)

        void store(String filename,
                  libcpp_vector[ProteinIdentification] & protein_ids,
                  libcpp_vector[PeptideIdentification] & peptide_ids,
                  )
