from libcpp cimport bool
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingTerm.h>" namespace "OpenMS":

    cdef cppclass CVMappingTerm:

        CVMappingTerm() except + nogil 
        CVMappingTerm(CVMappingTerm &) except + nogil 

        void setAccession(String accession) except + nogil  # wrap-doc:Sets the accession string of the term

        String getAccession() except + nogil  # wrap-doc:Returns the accession string of the term

        void setUseTermName(bool use_term_name) except + nogil  # wrap-doc:Sets whether the term name should be used, instead of the accession

        bool getUseTermName() except + nogil  # wrap-doc:Returns whether the term name should be used, instead of the accession

        void setUseTerm(bool use_term) except + nogil  # wrap-doc:Sets whether the term itself can be used (or only its children)

        bool getUseTerm() except + nogil  # wrap-doc:Returns true if the term can be used, false if only children are allowed

        void setTermName(String term_name) except + nogil  # wrap-doc:Sets the name of the term

        String getTermName() except + nogil  # wrap-doc:Returns the name of the term

        void setIsRepeatable(bool is_repeatable) except + nogil  # wrap-doc:Sets whether this term can be repeated

        bool getIsRepeatable() except + nogil  # wrap-doc:Returns true if this term can be repeated, false otherwise

        void setAllowChildren(bool allow_children) except + nogil  # wrap-doc:Sets whether children of this term are allowed

        bool getAllowChildren() except + nogil  # wrap-doc:Returns true if the children of this term are allowed to be used

        void setCVIdentifierRef(String cv_identifier_ref) except + nogil  # wrap-doc:Sets the CV identifier reference string, e.g. UO for unit obo

        String getCVIdentifierRef() except + nogil  # wrap-doc:Returns the CV identifier reference string

        # equality operator
        bool operator==(CVMappingTerm rhs) except + nogil 

        # inequality operator
        bool operator!=(CVMappingTerm rhs) except + nogil 
