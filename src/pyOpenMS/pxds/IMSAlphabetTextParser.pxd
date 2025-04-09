from Types cimport *
from IMSAlphabetParser cimport *

# cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>" namespace "OpenMS::ims":
#    
#    cdef cppclass IMSAlphabetTextParser(IMSAlphabetParser) :
#        # wrap-inherits:
#        #  IMSAlphabetParser
#        IMSAlphabetTextParser() except + nogil  
#        IMSAlphabetTextParser(IMSAlphabetTextParser &) except + nogil 

#        # ContainerType  getElements() except + nogil 
#        # NAMESPACE # void parse(std::istream & is_) except + nogil 
# 
