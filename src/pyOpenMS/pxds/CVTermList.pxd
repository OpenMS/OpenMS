from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from String cimport *
from MetaInfoInterface cimport *
from CVTerm cimport *
from CVMappingRule cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/METADATA/CVTermList.h>" namespace "OpenMS":

    cdef cppclass CVTermList(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        ######################################################################
        # Cython has a problem with inheritance of overloaded methods, so we
        # can only declare one of the two methods here.
        # see eg Precursor.pxd 

        CVTermList() except + nogil 
        CVTermList(CVTermList &) except + nogil 

        void setCVTerms(libcpp_vector[CVTerm] & terms)  except + nogil  # wrap-doc:Sets the CV terms
        void replaceCVTerm(CVTerm & term)               except + nogil  # wrap-doc:Replaces the specified CV term

        # will not wrap due to Cython inheritance issue
        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession) except + nogil 
        # void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map) except + nogil 

        void consumeCVTerms(libcpp_map[String, libcpp_vector[CVTerm] ] cv_term_map) except + nogil  # wrap-doc:Merges the given map into the member map, no duplicate checking

        libcpp_map[String, libcpp_vector[CVTerm] ] getCVTerms() except + nogil  # wrap-doc:Returns the accession string of the term
        void addCVTerm(CVTerm & term)                   except + nogil  # wrap-doc:Adds a CV term

        bool operator==(CVTermList)  except + nogil 
        bool operator!=(CVTermList)  except + nogil 

        bool hasCVTerm(String accession)  except + nogil 
        bool empty()                      except + nogil 

