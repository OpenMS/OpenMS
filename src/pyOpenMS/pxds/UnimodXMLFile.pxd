from Types cimport *
from XMLFile cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/FORMAT/UnimodXMLFile.h>" namespace "OpenMS":
    
    cdef cppclass UnimodXMLFile(XMLFile) :
        # wrap-inherits:
        #  XMLFile

        UnimodXMLFile() nogil except +
        UnimodXMLFile(UnimodXMLFile) nogil except + #wrap-ignore

        # TODO raw ptr in vector
        # void load(const String & filename, libcpp_vector[ ResidueModification * ] & modifications) nogil except +

