from libcpp cimport bool
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingTerm.h>" namespace "OpenMS":

    cdef cppclass CVMappingTerm:

        CVMappingTerm() nogil except +
        CVMappingTerm(CVMappingTerm &) nogil except +

        # sets the accession string of the term
        void setAccession(String accession) nogil except + # wrap-doc:Sets the accession string of the term

        # returns the accession string of the term
        String getAccession() nogil except + # wrap-doc:Returns the accession string of the term

        # sets whether the term name should be used, instead of the accession
        void setUseTermName(bool use_term_name) nogil except + # wrap-doc:Sets whether the term name should be used, instead of the accession

        # returns whether the term name should be used, instead of the accession
        bool getUseTermName() nogil except + # wrap-doc:Returns whether the term name should be used, instead of the accession

        # sets whether the term itself can be used (or only its children)
        void setUseTerm(bool use_term) nogil except + # wrap-doc:Sets whether the term itself can be used (or only its children)

        # returns true if the term can be used, false if only children are allowed
        bool getUseTerm() nogil except + # wrap-doc:Returns true if the term can be used, false if only children are allowed

        # sets the name of the term
        void setTermName(String term_name) nogil except + # wrap-doc:Sets the name of the term

        # returns the name of the term
        String getTermName() nogil except + # wrap-doc:Returns the name of the term

        # sets whether this term can be repeated
        void setIsRepeatable(bool is_repeatable) nogil except + # wrap-doc:Sets whether this term can be repeated

        # returns true if this term can be repeated, false otherwise
        bool getIsRepeatable() nogil except + # wrap-doc:Returns true if this term can be repeated, false otherwise

        # sets whether children of this term are allowed
        void setAllowChildren(bool allow_children) nogil except + # wrap-doc:Sets whether children of this term are allowed

        # returns true if the children of this term are allowed to be used
        bool getAllowChildren() nogil except + # wrap-doc:Returns true if the children of this term are allowed to be used

        # sets the cv identifier reference string, e.g. UO for unit obo
        void setCVIdentifierRef(String cv_identifier_ref) nogil except + # wrap-doc:Sets the CV identifier reference string, e.g. UO for unit obo

        # returns the cv identifier reference string
        String getCVIdentifierRef() nogil except + # wrap-doc:Returns the CV identifier reference string

        # equality operator
        bool operator==(CVMappingTerm rhs) nogil except +

        # inequality operator
        bool operator!=(CVMappingTerm rhs) nogil except +

