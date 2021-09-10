from libcpp cimport bool
from Map cimport *
from String cimport *
from MetaInfoInterface cimport *
from CVTerm cimport *
from CVMappingRule cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/METADATA/CVTermList.h>" namespace "OpenMS":

    cdef cppclass CVTermList(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        ######################################################################
        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the two methods here.
        # see eg Precursor.pxd 

        CVTermList() nogil except +
        CVTermList(CVTermList &) nogil except +

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except + # wrap-doc:Sets the CV terms
        void replaceCVTerm(CVTerm & term)               nogil except + # wrap-doc:Replaces the specified CV term

        # will not wrap due to Cython inheritance issue
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        # void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +

        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except + # wrap-doc:Merges the given map into the member map, no duplicate checking

        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except + # wrap-doc:Returns the accession string of the term
        void addCVTerm(CVTerm & term)                   nogil except + # wrap-doc:Adds a CV term

        bool operator==(CVTermList)  nogil except +
        bool operator!=(CVTermList)  nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

