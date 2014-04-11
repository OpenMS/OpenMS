from Types cimport *
from libcpp cimport bool
from String cimport *

cdef extern from "<OpenMS/FORMAT/XMLFile.h>" namespace "OpenMS::Internal":
    
    cdef cppclass XMLFile "OpenMS::Internal::XMLFile":
        XMLFile() nogil except +
        XMLFile(XMLFile) nogil except + #wrap-ignore
        XMLFile(String & schema_location, String & version) nogil except +
        # NAMESPACE # bool isValid(String & filename, std::ostream & os) nogil except +
        String  getVersion() nogil except +

