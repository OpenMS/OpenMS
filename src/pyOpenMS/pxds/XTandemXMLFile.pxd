from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from ModificationDefinitionsSet cimport *

cdef extern from "<OpenMS/FORMAT/XTandemXMLFile.h>" namespace "OpenMS":

    cdef cppclass XTandemXMLFile:

        XTandemXMLFile() nogil except +

        void load(String filename, ProteinIdentification & protein_identification, libcpp_vector[PeptideIdentification] & id_data) nogil except +

        # /// sets the valid modifications
        void setModificationDefinitionsSet(ModificationDefinitionsSet rhs) nogil except +

