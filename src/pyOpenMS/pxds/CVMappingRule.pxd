from libcpp cimport bool
from String cimport *
from CVTerm cimport *
from CVMappingTerm cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingRule.h>" namespace "OpenMS":

    cdef cppclass CVMappingRule:

        CVMappingRule() nogil except +
        CVMappingRule(CVMappingRule &) nogil except +

        void setIdentifier(String identifier) nogil except + # wrap-doc:Sets the identifier of the rule

        String getIdentifier() nogil except + # wrap-doc:Returns the identifier of the rule

        void setElementPath(String element_path) nogil except + # wrap-doc:Sets the path of the DOM element, where this rule is allowed

        String getElementPath() nogil except + # wrap-doc:Returns the path of the DOM element, where this rule is allowed

        void setRequirementLevel(RequirementLevel level) nogil except + # wrap-doc:Sets the requirement level of this rule

        RequirementLevel getRequirementLevel() nogil except + # wrap-doc:Returns the requirement level of this rule

        void setCombinationsLogic(CombinationsLogic combinations_logic) nogil except + # wrap-doc:Sets the combination operator of the rule

        CombinationsLogic getCombinationsLogic() nogil except + # wrap-doc:Returns the combinations operator of the rule

        void setScopePath(String path) nogil except + # wrap-doc:Sets the scope path of the rule

        String getScopePath() nogil except + # wrap-doc:Returns the scope path of the rule

        void setCVTerms(libcpp_vector[CVMappingTerm] cv_terms) nogil except + # wrap-doc:Sets the terms which are allowed

        libcpp_vector[CVMappingTerm] getCVTerms() nogil except + # wrap-doc:Returns the allowed terms

        void addCVTerm(CVMappingTerm cv_terms) nogil except + # wrap-doc:Adds a term to the allowed terms

        # equality operator
        bool operator==(CVMappingRule rhs) nogil except +

        # inequality operator
        bool operator!=(CVMappingRule rhs) nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingRule.h>" namespace "OpenMS::CVMappingRule":

    # enum to specify the requirement level
    cdef enum RequirementLevel:
        # wrap-attach:
        #     CVMappingRule
        MUST = 0,
        SHOULD = 1,
        MAY = 2

    # enum to specify the combination operator
    cdef enum CombinationsLogic:
        # wrap-attach:
        #     CVMappingRule
        OR = 0,
        AND = 1,
        XOR = 2
