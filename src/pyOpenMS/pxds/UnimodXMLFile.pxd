from Types cimport *
from XMLFile cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/FORMAT/UnimodXMLFile.h>" namespace "OpenMS":
    
    cdef cppclass UnimodXMLFile(XMLFile) :
        # wrap-doc:
        #  Used to load XML files from unimod.org files
        # wrap-inherits:
        #  XMLFile

        UnimodXMLFile() except + nogil 
        # private
        UnimodXMLFile(UnimodXMLFile) except + nogil  # wrap-ignore

        # TODO raw ptr in vector
        # void load(const String & filename, libcpp_vector[ ResidueModification * ] & modifications) except + nogil 
