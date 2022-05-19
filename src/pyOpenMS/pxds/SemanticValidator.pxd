from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from StringList cimport *
from CVMappings cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>" namespace "OpenMS::Internal":

    cdef cppclass SemanticValidator:

        # private
        SemanticValidator() nogil except + # wrap-ignore
        # private
        SemanticValidator(SemanticValidator &) nogil except +  # wrap-ignore
        SemanticValidator(CVMappings mapping, ControlledVocabulary cv) nogil except +

        bool validate(String filename, StringList errors, StringList warnings) nogil except +

        bool locateTerm(String path, SemanticValidator_CVTerm & parsed_term) nogil except + # wrap-doc:Checks if a CVTerm is allowed in a given path

        void setTag(String tag) nogil except + # wrap-doc:Sets the CV parameter tag name (default 'cvParam')

        void setAccessionAttribute(String accession) nogil except + # wrap-doc:Sets the name of the attribute for accessions in the CV parameter tag name (default 'accession')

        void setNameAttribute(String name) nogil except + # wrap-doc:Sets the name of the attribute for accessions in the CV parameter tag name (default 'name')

        void setValueAttribute(String value) nogil except + # wrap-doc:Sets the name of the attribute for accessions in the CV parameter tag name (default 'value')

        void setCheckTermValueTypes(bool check) nogil except + # wrap-doc:Sets if CV term value types should be check (enabled by default)

        void setCheckUnits(bool check) nogil except + # wrap-doc:Sets if CV term units should be check (disabled by default)

        void setUnitAccessionAttribute(String accession) nogil except + # wrap-doc:Sets the name of the unit accession attribute (default 'unitAccession')

        void setUnitNameAttribute(String name) nogil except + # wrap-doc:Sets the name of the unit name attribute (default 'unitName')

cdef extern from "<OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>" namespace "OpenMS::Internal::SemanticValidator":

      cdef cppclass SemanticValidator_CVTerm "OpenMS::Internal::SemanticValidator::CVTerm":
        SemanticValidator_CVTerm() nogil except + # compiler
        SemanticValidator_CVTerm(SemanticValidator_CVTerm &) nogil except + # compiler
        String accession
        String name
        String value
        bool has_value
        String unit_accession
        bool has_unit_accession
        String unit_name
        bool has_unit_name
