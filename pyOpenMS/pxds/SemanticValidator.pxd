from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from StringList cimport *
from CVMappings cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>" namespace "OpenMS::Internal":

    cdef cppclass SemanticValidator:

        SemanticValidator(CVMappings mapping, ControlledVocabulary cv) nogil except +

        bool validate(String filename, StringList errors, StringList warnings) nogil except +

        # Checks if a CVTerm is allowed in a given path
        # TODO yet anothre CV Term
        # bool locateTerm(String path, SemanticValidator_CVTerm & parsed_term) nogil except +

        # Sets the CV parameter tag name (default: 'cvParam')
        void setTag(String tag) nogil except +

        # Sets the name of the attribute for accessions in the CV parameter tag name (default: 'accession')
        void setAccessionAttribute(String accession) nogil except +

        # Sets the name of the attribute for accessions in the CV parameter tag name (default: 'name')
        void setNameAttribute(String name) nogil except +

        # Sets the name of the attribute for accessions in the CV parameter tag name (default: 'value')
        void setValueAttribute(String value) nogil except +

        # Set if CV term value types should be check (enabled by default)
        void setCheckTermValueTypes(bool check) nogil except +

        # Set if CV term units should be check (disabled by default)
        void setCheckUnits(bool check) nogil except +

        # Sets the name of the unit accession attribute (default: 'unitAccession')
        void setUnitAccessionAttribute(String accession) nogil except +

        # Sets the name of the unit name attribute (default: 'unitName')
        void setUnitNameAttribute(String name) nogil except +

