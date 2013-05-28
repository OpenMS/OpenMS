from Types cimport *
from libcpp.string cimport string as libcpp_string

# cannot wrap template argument std::istream
# typename Container = std::map<std::string, double>,
# cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>" namespace "OpenMS::ims":
#     
#     cdef cppclass IMSAlphabetParser: # [AlphabetElementType,Container,InputSource]:
#         IMSAlphabetParser(IMSAlphabetParser) nogil except + #wrap-ignore
#         void load(libcpp_string & fname) nogil except +
#         # ContainerType getElements() nogil except +
#         # void parse(InputSource & is_) nogil except +

