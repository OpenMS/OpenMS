from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from MetaInfoInterface cimport *
from CVTerm cimport *
from Map cimport *

cdef extern from "<OpenMS/METADATA/CVTermListInterface.h>" namespace "OpenMS":
    
    cdef cppclass CVTermListInterface(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface

        CVTermListInterface() nogil except +
        CVTermListInterface(CVTermListInterface) nogil except +

        bool operator==(CVTermListInterface & rhs) nogil except +
        bool operator!=(CVTermListInterface & rhs) nogil except +

        void replaceCVTerms(Map[ String, libcpp_vector[ CVTerm ] ] & cv_terms) nogil except +
        void setCVTerms(libcpp_vector[ CVTerm ] & terms) nogil except +

        ######################################################################
        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the two methods here.
        void replaceCVTerm(CVTerm & cv_term) nogil except +
        void replaceCVTerms(libcpp_vector[ CVTerm ] & cv_terms, const String & accession) nogil except +
        # void replaceCVTerms(Map[ String, libcpp_vector[ CVTerm ] ] & cv_term_map) nogil except +

        void consumeCVTerms(Map[ String, libcpp_vector[ CVTerm ] ] & cv_term_map) nogil except +
        Map[ String, libcpp_vector[ CVTerm ] ]  getCVTerms() nogil except +


        void addCVTerm(CVTerm & term) nogil except +
        bool hasCVTerm(const String & accession) nogil except +
        bool empty() nogil except +

