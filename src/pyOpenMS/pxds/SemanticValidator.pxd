from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from StringList cimport *
from CVMappings cimport *
from ControlledVocabulary cimport *

cdef extern from "<OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>" namespace "OpenMS::Internal":

    cdef cppclass SemanticValidator:

        # private
        SemanticValidator() except + nogil  # wrap-ignore
        # private
        SemanticValidator(SemanticValidator &) except + nogil   # wrap-ignore
        SemanticValidator(CVMappings mapping, ControlledVocabulary cv) except + nogil 

        bool validate(String filename, StringList errors, StringList warnings) except + nogil 

        bool locateTerm(String path, SemanticValidator_CVTerm & parsed_term) except + nogil  # wrap-doc:Checks if a CVTerm is allowed in a given path

        void setTag(String tag) except + nogil  # wrap-doc:Sets the CV parameter tag name (default 'cvParam')

        void setAccessionAttribute(String accession) except + nogil  # wrap-doc:Sets the name of the attribute for accessions in the CV parameter tag name (default 'accession')

        void setNameAttribute(String name) except + nogil  # wrap-doc:Sets the name of the attribute for accessions in the CV parameter tag name (default 'name')

        void setValueAttribute(String value) except + nogil  # wrap-doc:Sets the name of the attribute for accessions in the CV parameter tag name (default 'value')

        void setCheckTermValueTypes(bool check) except + nogil  # wrap-doc:Sets if CV term value types should be check (enabled by default)

        void setCheckUnits(bool check) except + nogil  # wrap-doc:Sets if CV term units should be check (disabled by default)

        void setUnitAccessionAttribute(String accession) except + nogil  # wrap-doc:Sets the name of the unit accession attribute (default 'unitAccession')

        void setUnitNameAttribute(String name) except + nogil  # wrap-doc:Sets the name of the unit name attribute (default 'unitName')

cdef extern from "<OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>" namespace "OpenMS::Internal::SemanticValidator":

      cdef cppclass SemanticValidator_CVTerm "OpenMS::Internal::SemanticValidator::CVTerm":
        SemanticValidator_CVTerm() except + nogil  # compiler
        SemanticValidator_CVTerm(SemanticValidator_CVTerm &) except + nogil  # compiler
        String accession
        String name
        String value
        bool has_value
        String unit_accession
        bool has_unit_accession
        String unit_name
        bool has_unit_name
