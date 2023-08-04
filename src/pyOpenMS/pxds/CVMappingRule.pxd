from libcpp cimport bool
from String cimport *
from CVTerm cimport *
from CVMappingTerm cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingRule.h>" namespace "OpenMS":

    cdef cppclass CVMappingRule:

        CVMappingRule() except + nogil 
        CVMappingRule(CVMappingRule &) except + nogil 

        void setIdentifier(String identifier) except + nogil  # wrap-doc:Sets the identifier of the rule

        String getIdentifier() except + nogil  # wrap-doc:Returns the identifier of the rule

        void setElementPath(String element_path) except + nogil  # wrap-doc:Sets the path of the DOM element, where this rule is allowed

        String getElementPath() except + nogil  # wrap-doc:Returns the path of the DOM element, where this rule is allowed

        void setRequirementLevel(RequirementLevel level) except + nogil  # wrap-doc:Sets the requirement level of this rule

        RequirementLevel getRequirementLevel() except + nogil  # wrap-doc:Returns the requirement level of this rule

        void setCombinationsLogic(CombinationsLogic combinations_logic) except + nogil  # wrap-doc:Sets the combination operator of the rule

        CombinationsLogic getCombinationsLogic() except + nogil  # wrap-doc:Returns the combinations operator of the rule

        void setScopePath(String path) except + nogil  # wrap-doc:Sets the scope path of the rule

        String getScopePath() except + nogil  # wrap-doc:Returns the scope path of the rule

        void setCVTerms(libcpp_vector[CVMappingTerm] cv_terms) except + nogil  # wrap-doc:Sets the terms which are allowed

        libcpp_vector[CVMappingTerm] getCVTerms() except + nogil  # wrap-doc:Returns the allowed terms

        void addCVTerm(CVMappingTerm cv_terms) except + nogil  # wrap-doc:Adds a term to the allowed terms

        # equality operator
        bool operator==(CVMappingRule rhs) except + nogil 

        # inequality operator
        bool operator!=(CVMappingRule rhs) except + nogil 

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappingRule.h>" namespace "OpenMS::CVMappingRule":

    # enum to specify the requirement level
    cdef enum RequirementLevel:
        # wrap-attach:
        #    CVMappingRule
        MUST = 0,
        SHOULD = 1,
        MAY = 2

    # enum to specify the combination operator
    cdef enum CombinationsLogic:
        # wrap-attach:
        #    CVMappingRule
        OR = 0,
        AND = 1,
        XOR = 2
