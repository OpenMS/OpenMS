from Types cimport *
from IMSAlphabetParser cimport *

# cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>" namespace "OpenMS::ims":
#     
#     cdef cppclass IMSAlphabetTextParser(IMSAlphabetParser) :
#         # wrap-inherits:
#         #  IMSAlphabetParser
#         IMSAlphabetTextParser() nogil except + 
#         IMSAlphabetTextParser(IMSAlphabetTextParser) nogil except + #wrap-ignore
#         # ContainerType  getElements() nogil except +
#         # NAMESPACE # void parse(std::istream & is_) nogil except +
# 
