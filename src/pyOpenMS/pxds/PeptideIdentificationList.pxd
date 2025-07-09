from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *

cdef extern from "<OpenMS/METADATA/PeptideIdentificationList.h>" namespace "OpenMS":

    cdef cppclass PeptideIdentificationList:
        # wrap-doc:
        #  Container for peptide identifications from multiple spectra.

        PeptideIdentificationList() except + nogil 
        PeptideIdentificationList(PeptideIdentificationList) except + nogil 
        
        # Vector-like interface
        int size() except + nogil 
        bool empty() except + nogil 
        void clear() except + nogil 
        void push_back(PeptideIdentification) except + nogil 
        
        # Element access
        PeptideIdentification operator[](size_t) except + nogil 
        PeptideIdentification at(size_t) except + nogil 
        PeptideIdentification back() except + nogil 
        PeptideIdentification front() except + nogil 
        
        # Iterators for Python iteration support
        libcpp_vector[PeptideIdentification].iterator begin() except + nogil     # wrap-iter-begin:__iter__(PeptideIdentification)
        libcpp_vector[PeptideIdentification].iterator end()   except + nogil     # wrap-iter-end:__iter__(PeptideIdentification)