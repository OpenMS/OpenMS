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
        # cython has a problem with inheritance of overloaded methods,
        # see eg Precursor.pxd 

        CVTermList()            nogil except +
        CVTermList(CVTermList)  nogil except +

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) nogil except +
        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +

        void consumeCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except +
        void addCVTerm(CVTerm & term)                   nogil except +

        bool operator==(CVTermList)  nogil except +
        bool operator!=(CVTermList)  nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +

