from CVTerm cimport *
from Map cimport *
from String cimport *
from MetaInfoInterface cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/METADATA/CVTermList.h>" namespace "OpenMS":

    cdef cppclass CVTermList(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        CVTermList()            nogil except +
        CVTermList(CVTermList)  nogil except +

        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) nogil except +

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool operator==(CVTermList)  nogil except +
        bool operator!=(CVTermList)  nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +


