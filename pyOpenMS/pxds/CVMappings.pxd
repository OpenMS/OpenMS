from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from CVReference cimport *
from CVMappingRule cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappings.h>" namespace "OpenMS":

    cdef cppclass CVMappings:

        CVMappings() nogil except +
        CVMappings(CVMappings) nogil except +
        void setMappingRules(libcpp_vector[ CVMappingRule ] &cv_mapping_rules)
        libcpp_vector[ CVMappingRule ]  getMappingRules()
        void addMappingRule(CVMappingRule &cv_mapping_rule)
        void setCVReferences(libcpp_vector[ CVReference ] &cv_references)
        libcpp_vector[ CVReference ]  getCVReferences()
        void addCVReference(CVReference &cv_reference)
        bool hasCVReference(String &identifier)

