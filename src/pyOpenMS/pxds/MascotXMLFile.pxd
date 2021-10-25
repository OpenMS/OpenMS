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
        MascotXMLFile(MascotXMLFile &) nogil except +

        void load(const String & filename,
                  ProteinIdentification & protein_identification,
                  libcpp_vector[ PeptideIdentification ] & id_data,
                  SpectrumMetaDataLookup & rt_mapping) nogil except +
            # wrap-doc:
                #   Loads data from a Mascot XML file
                #   -----
                #   :param filename: The file to be loaded
                #   :param protein_identification: Protein identifications belonging to the whole experiment
                #   :param id_data: The identifications with m/z and RT
                #   :param lookup: Helper object for looking up spectrum meta data
                #   :raises:
                #     Exception: FileNotFound is thrown if the file does not exists
                #   :raises:
                #     Exception: ParseError is thrown if the file does not suit to the standard

        # TODO fix
        # void load(const String & filename,
        #           ProteinIdentification & protein_identification,
        #           libcpp_vector[ PeptideIdentification ] & id_data,
        #           libcpp_map[ String, libcpp_vector[ AASequence ] ] & peptides,
        #           SpectrumMetaDataLookup & rt_mapping) nogil except +

        void initializeLookup(SpectrumMetaDataLookup & lookup, MSExperiment& experiment, const String & scan_regex) nogil except +
            # wrap-doc:
                #   Initializes a helper object for looking up spectrum meta data (RT, m/z)
                #   -----
                #   :param lookup: Helper object to initialize
                #   :param experiment: Experiment containing the spectra
                #   :param scan_regex: Optional regular expression for extracting information from references to spectra

