from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *

cdef extern from "<OpenMS/METADATA/PeptideIdentificationList.h>" namespace "OpenMS":

    cdef cppclass PeptideIdentificationList(libcpp_vector[PeptideIdentification]):
        # wrap-inherits:
        #   libcpp_vector[PeptideIdentification]
        # wrap-doc:
        #  Container for peptide identifications from multiple spectra
        #  
        #  This class represents a collection of PeptideIdentification objects,
        #  typically from multiple spectra in a single identification run.
        #  It provides all the functionality of std::vector<PeptideIdentification>
        #  while maintaining type safety and allowing for future extensions.

        PeptideIdentificationList() except + nogil 
        PeptideIdentificationList(PeptideIdentificationList &) except + nogil 
        PeptideIdentificationList(libcpp_vector[PeptideIdentification] &) except + nogil
        
        # Note: Most vector functionality is automatically wrapped through inheritance
        # Additional methods can be added here if needed