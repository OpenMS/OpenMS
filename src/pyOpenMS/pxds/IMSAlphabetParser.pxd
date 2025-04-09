from Types cimport *
from libcpp.string cimport string as libcpp_string

# cannot wrap template argument std::istream
# typename Container = std::map<std::string, double>,
# cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>" namespace "OpenMS::ims":
#    
#    cdef cppclass IMSAlphabetParser: # [AlphabetElementType,Container,InputSource]:
#        IMSAlphabetParser() except + nogil  # compiler
#        IMSAlphabetParser(IMSAlphabetParser &) except + nogil  # compiler
#        void load(libcpp_string & fname) except + nogil 
#        # ContainerType getElements() except + nogil 
#        # void parse(InputSource & is_) except + nogil 
