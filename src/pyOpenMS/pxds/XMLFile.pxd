from Types cimport *
from libcpp cimport bool
from String cimport *

cdef extern from "<OpenMS/FORMAT/XMLFile.h>" namespace "OpenMS::Internal":
    
    cdef cppclass XMLFile "OpenMS::Internal::XMLFile":
        XMLFile() except + nogil 
        XMLFile(XMLFile &) except + nogil  

        XMLFile(const String & schema_location, const String & version) except + nogil 
        # NAMESPACE # bool isValid(const String & filename, std::ostream & os) except + nogil 
        String getVersion() except + nogil  # wrap-doc:Return the version of the schema
