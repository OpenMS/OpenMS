from CVTerm cimport *
from Map cimport *
from String cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/CVTermList.h>" namespace "OpenMS":

    cdef cppclass CVTermList(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface
         CVTermList()
         # fails with ! 'raise Exception("can not handle wrapped classes as keys in map")'
         # Map[String, libcpp_vector[CVTerm] ] getCVTerms()
         void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
         void addCVTerm(CVTerm & term)                   nogil except +


