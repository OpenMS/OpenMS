from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *
from CVTerm cimport *


cdef extern from "<OpenMS/METADATA/CVTermListInterface.h>" namespace "OpenMS":
    
    cdef cppclass CVTermListInterface(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface

        CVTermListInterface() except + nogil 
        CVTermListInterface(CVTermListInterface &) except + nogil 

        bool operator==(CVTermListInterface & rhs) except + nogil 
        bool operator!=(CVTermListInterface & rhs) except + nogil 

        void replaceCVTerms(libcpp_map[ String, libcpp_vector[ CVTerm ] ] & cv_terms) except + nogil 
        void setCVTerms(libcpp_vector[ CVTerm ] & terms) except + nogil 

        ######################################################################
        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the two methods here.
        void replaceCVTerm(CVTerm & cv_term) except + nogil 
        void replaceCVTerms(libcpp_vector[ CVTerm ] & cv_terms, const String & accession) except + nogil 
        # void replaceCVTerms(Map[ String, libcpp_vector[ CVTerm ] ] & cv_term_map) except + nogil 

        void consumeCVTerms(libcpp_map[ String, libcpp_vector[ CVTerm ] ] & cv_term_map) except + nogil  # wrap-doc:Merges the given map into the member map, no duplicate checking
        libcpp_map[ String, libcpp_vector[ CVTerm ] ]  getCVTerms() except + nogil 


        void addCVTerm(CVTerm & term) except + nogil  # wrap-doc:Adds a CV term
        bool hasCVTerm(const String & accession) except + nogil  # wrap-doc:Checks whether the term has a value
        bool empty() except + nogil 

