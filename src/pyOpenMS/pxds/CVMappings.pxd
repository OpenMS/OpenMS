from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from CVReference cimport *
from CVMappingRule cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappings.h>" namespace "OpenMS":

    cdef cppclass CVMappings:

        CVMappings() except + nogil 
        CVMappings(CVMappings &) except + nogil 
        void setMappingRules(libcpp_vector[ CVMappingRule ]& cv_mapping_rules) except + nogil  # wrap-doc:Sets the mapping rules of the mapping file
        libcpp_vector[ CVMappingRule ]  getMappingRules() except + nogil  # wrap-doc:Returns the mapping rules
        void addMappingRule(CVMappingRule& cv_mapping_rule) except + nogil  # wrap-doc:Adds a mapping rule
        void setCVReferences(libcpp_vector[ CVReference ]& cv_references) except + nogil  # wrap-doc:Sets the CV references
        libcpp_vector[ CVReference ] getCVReferences() except + nogil  # wrap-doc:Returns the CV references
        void addCVReference(CVReference& cv_reference) except + nogil  # wrap-doc:Adds a CV reference
        bool hasCVReference(const String& identifier) except + nogil  # wrap-doc:Returns true if a CV reference is given
