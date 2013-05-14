from libcpp cimport bool
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingTerm.h>" namespace "OpenMS":

    cdef cppclass CVMappingTerm:

        CVMappingTerm()               nogil except +
        CVMappingTerm(CVMappingTerm)  nogil except +

        # sets the accession string of the term
        void setAccession(String accession) nogil except +

        # returns the accession string of the term
        String getAccession() nogil except +

        # sets whether the term name should be used, instead of the accession
        void setUseTermName(bool use_term_name) nogil except +

        # returns whether the term name should be used, instead of the accession
        bool getUseTermName() nogil except +

        # sets whether the term itself can be used (or only its children)
        void setUseTerm(bool use_term) nogil except +

        # returns true if the term can be used, false if only children are allowed
        bool getUseTerm() nogil except +

        # sets the name of the term
        void setTermName(String term_name) nogil except +

        # returns the name of the term
        String getTermName() nogil except +

        # sets whether this term can be repeated
        void setIsRepeatable(bool is_repeatable) nogil except +

        # returns true if this term can be repeated, false otherwise
        bool getIsRepeatable() nogil except +

        # sets whether children of this term are allowed
        void setAllowChildren(bool allow_children) nogil except +

        # returns true if the children of this term are allowed to be used
        bool getAllowChildren() nogil except +

        # sets the cv identifier reference string, e.g. UO for unit obo
        void setCVIdentifierRef(String cv_identifier_ref) nogil except +

        # returns the cv identifier reference string
        String getCVIdentifierRef() nogil except +

        # equality operator
        bool operator==(CVMappingTerm rhs) nogil except +

        # inequality operator
        bool operator!=(CVMappingTerm rhs) nogil except +

