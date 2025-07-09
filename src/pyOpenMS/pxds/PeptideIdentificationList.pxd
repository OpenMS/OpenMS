from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *

cdef extern from "<OpenMS/METADATA/PeptideIdentificationList.h>" namespace "OpenMS":

    cdef cppclass PeptideIdentificationList:
        # wrap-doc:
        #  Container for peptide identifications from multiple spectra.
        #  
        #  This class represents a collection of PeptideIdentification objects,
        #  typically from multiple spectra in a single identification run.
        #  It provides all the functionality of std::vector<PeptideIdentification>
        #  while maintaining type safety and allowing for future extensions.
        #  
        #  Uses composition over inheritance to avoid the pitfalls of inheriting
        #  from STL containers.

        PeptideIdentificationList() except + nogil 
        PeptideIdentificationList(PeptideIdentificationList &) except + nogil 
        PeptideIdentificationList(libcpp_vector[PeptideIdentification] &) except + nogil
        PeptideIdentificationList(size_t) except + nogil 
        PeptideIdentificationList(size_t, PeptideIdentification &) except + nogil 
        
        # Vector-like interface
        int size() except + nogil 
        bool empty() except + nogil 
        void clear() except + nogil 
        void push_back(PeptideIdentification &) except + nogil 
        
        # Element access
        PeptideIdentification & operator[](size_t) except + nogil 
        PeptideIdentification & at(size_t) except + nogil 
        PeptideIdentification & back() except + nogil 
        PeptideIdentification & front() except + nogil 
        
        # Iterators for Python iteration support
        libcpp_vector[PeptideIdentification].iterator begin() except + nogil     # wrap-iter-begin:__iter__(PeptideIdentification)
        libcpp_vector[PeptideIdentification].iterator end()   except + nogil     # wrap-iter-end:__iter__(PeptideIdentification)
        
        # Conversion operators for compatibility
        # Note: These are handled automatically by Cython when needed