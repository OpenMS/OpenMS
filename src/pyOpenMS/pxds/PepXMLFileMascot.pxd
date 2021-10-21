from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from XMLHandler cimport *
from XMLFile cimport *
from AASequence cimport *

cdef extern from "<OpenMS/FORMAT/PepXMLFileMascot.h>" namespace "OpenMS":
    
    cdef cppclass PepXMLFileMascot :
        # wrap-doc:
            #   Used to load Mascot PepXML files
            #   -----
            #   A schema for this format can be found at http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd

        PepXMLFileMascot() nogil except +
        # copy constructor of 'PepXMLFileMascot' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor protected Internal::XMLHandler,
        PepXMLFileMascot(PepXMLFileMascot &) nogil except + # wrap-ignore

        # TODO map
        # void load(const String & filename, libcpp_map[ String, libcpp_vector[ AASequence ] ] & peptides) nogil except +

