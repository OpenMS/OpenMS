from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from PeptideIdentification cimport *

cdef extern from "<OpenMS/METADATA/PeptideIdentificationList.h>" namespace "OpenMS":

    cdef cppclass PeptideIdentificationList:
        # wrap-doc:
        #  A container for peptide identifications from multiple spectra.
        #
        #  This class provides a vector-like interface for managing collections
        #  of peptide identifications, typically obtained from tandem mass
        #  spectrometry experiments or peptide database searches. Each PeptideIdentification represents
        #  the identification results from a single spectrum.
        #
        #  This class supports direct iteration in Python.

        PeptideIdentificationList() except + nogil 
        PeptideIdentificationList(PeptideIdentificationList) except + nogil 
        
        # Vector-like interface
        int size() except + nogil  # wrap-doc:Returns the number of peptide identifications
        bool empty() except + nogil  # wrap-doc:Returns true if the container is empty
        void clear() except + nogil  # wrap-doc:Removes all peptide identifications from the container
        void push_back(PeptideIdentification) except + nogil  # wrap-doc:Adds a peptide identification to the end of the container
        
        # Element access
        PeptideIdentification & operator[](size_t)      except + nogil  #wrap-upper-limit:size()
        PeptideIdentification at(size_t) except + nogil  # wrap-doc:Returns the peptide identification at the given index with bounds checking
        PeptideIdentification back() except + nogil  # wrap-doc:Returns the last peptide identification in the container
        PeptideIdentification front() except + nogil  # wrap-doc:Returns the first peptide identification in the container
        
        # Iterators for Python iteration support
        libcpp_vector[PeptideIdentification].iterator begin() except + nogil     # wrap-iter-begin:__iter__(PeptideIdentification)
        libcpp_vector[PeptideIdentification].iterator end()   except + nogil     # wrap-iter-end:__iter__(PeptideIdentification)