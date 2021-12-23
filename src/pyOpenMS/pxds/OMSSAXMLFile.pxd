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
        # wrap-doc:
        #   Used to load OMSSAXML files
        #   -----
        #   This class is used to load documents that implement
        #   the schema of OMSSAXML files

        OMSSAXMLFile() nogil except +
        # private
        OMSSAXMLFile(OMSSAXMLFile &) nogil except + # wrap-ignore

        void load(const String & filename,
                  ProteinIdentification & protein_identification,
                  libcpp_vector[ PeptideIdentification ] & id_data,
                  bool load_proteins,
                  bool load_empty_hits) nogil except +
            # wrap-doc:
                #   Loads data from a OMSSAXML file
                #   -----
                #   :param filename: The file to be loaded
                #   :param protein_identification: Protein identifications belonging to the whole experiment
                #   :param id_data: The identifications with m/z and RT
                #   :param load_proteins: If this flag is set to false, the protein identifications are not loaded
                #   :param load_empty_hits: Many spectra will not return a hit. Report empty peptide identifications?

        void setModificationDefinitionsSet(ModificationDefinitionsSet rhs) nogil except + # wrap-doc:Sets the valid modifications


