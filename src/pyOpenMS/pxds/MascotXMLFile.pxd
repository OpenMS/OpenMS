from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from AASequence cimport *
from XMLFile cimport *
from SpectrumMetaDataLookup cimport *
from Map cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/MascotXMLFile.h>" namespace "OpenMS":
    
    cdef cppclass MascotXMLFile(XMLFile) :
        # wrap-inherits:
        #  XMLFile
        MascotXMLFile() nogil except +
        MascotXMLFile(MascotXMLFile) nogil except + #wrap-ignore

        void load(const String & filename,
                  ProteinIdentification & protein_identification,
                  libcpp_vector[ PeptideIdentification ] & id_data,
                  SpectrumMetaDataLookup & rt_mapping) nogil except +

        # TODO fix
        # void load(const String & filename,
        #           ProteinIdentification & protein_identification,
        #           libcpp_vector[ PeptideIdentification ] & id_data,
        #           libcpp_map[ String, libcpp_vector[ AASequence ] ] & peptides,
        #           SpectrumMetaDataLookup & rt_mapping) nogil except +

        void initializeLookup(SpectrumMetaDataLookup & lookup, MSExperiment& experiment, const String & scan_regex) nogil except +

