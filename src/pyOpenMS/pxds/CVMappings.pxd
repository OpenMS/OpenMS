from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from CVReference cimport *
from CVMappingRule cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/CVMappings.h>" namespace "OpenMS":

    cdef cppclass CVMappings:

        CVMappings() nogil except +
        CVMappings(CVMappings) nogil except +
        void setMappingRules(libcpp_vector[ CVMappingRule ]& cv_mapping_rules) nogil except +
        libcpp_vector[ CVMappingRule ]  getMappingRules() nogil except +
        void addMappingRule(CVMappingRule& cv_mapping_rule) nogil except +
        void setCVReferences(libcpp_vector[ CVReference ]& cv_references) nogil except +
        libcpp_vector[ CVReference ] getCVReferences() nogil except +
        void addCVReference(CVReference& cv_reference) nogil except +
        bool hasCVReference(const String& identifier) nogil except +

