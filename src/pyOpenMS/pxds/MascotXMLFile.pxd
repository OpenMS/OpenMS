from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *
from AASequence cimport *
from XMLFile cimport *
from Map cimport *
from MSExperiment cimport *

# typedef Map<Size, float> RTMapping;
cdef extern from "<OpenMS/FORMAT/MascotXMLFile.h>" namespace "OpenMS":
    
    cdef cppclass MascotXMLFile(XMLFile) :
        # wrap-inherits:
        #  XMLFile
        MascotXMLFile() nogil except +
        MascotXMLFile(MascotXMLFile) nogil except + #wrap-ignore

        # TODO get RT mapping to work
        # void load(String & filename,
        #           ProteinIdentification & protein_identification,
        #           libcpp_vector[ PeptideIdentification ] & id_data,
        #           Map[size_t,float] & rt_mapping, String & scan_regex) nogil except +

        # void load(String & filename,
        #           ProteinIdentification & protein_identification,
        #           libcpp_vector[ PeptideIdentification ] & id_data,
        #           libcpp_map[ String, libcpp_vector[ AASequence ] ] & peptides,
        #           RTMapping & rt_mapping, String & scan_regex) nogil except +

        # NAMESPACE # void generateRTMapping(MSExperiment[Peak1D, ChromatogramPeak]::ConstIterator begin, MSExperiment[Peak1D, ChromatogramPeak]::ConstIterator end, RTMapping & rt_mapping) nogil except +

