from Types cimport *
from XMLFile cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from ResidueModification cimport *
from ModificationDefinitionsSet cimport *

cdef extern from "<OpenMS/FORMAT/OMSSAXMLFile.h>" namespace "OpenMS":
    
    cdef cppclass OMSSAXMLFile(XMLFile) :
        # wrap-inherits:
        #  XMLFile

        OMSSAXMLFile() nogil except +
        OMSSAXMLFile(OMSSAXMLFile) nogil except + #wrap-ignore

        void load(const String & filename,
                  ProteinIdentification & protein_identification,
                  libcpp_vector[ PeptideIdentification ] & id_data,
                  bool load_proteins,
                  bool load_empty_hits) nogil except +

        void setModificationDefinitionsSet(ModificationDefinitionsSet rhs) nogil except +


