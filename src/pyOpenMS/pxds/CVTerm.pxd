from libcpp cimport bool
from Types cimport *
from String cimport *
from DataValue cimport *

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS":

    cdef cppclass CVTerm:
         CVTerm() except + nogil 
         CVTerm(CVTerm &) except + nogil 

         bool operator==(CVTerm) except + nogil 

         void setAccession(String accession) except + nogil  # wrap-doc:Sets the accession string of the term
         String getAccession() except + nogil  # wrap-doc:Returns the accession string of the term

         void setName(String name) except + nogil  # wrap-doc:Sets the name of the term
         String getName() except + nogil  # wrap-doc:Returns the name of the term

         void setCVIdentifierRef(String cv_id_ref) except + nogil  # wrap-doc:Sets the CV identifier reference string, e.g. UO for unit obo
         String getCVIdentifierRef() except + nogil  # wrap-doc:Returns the CV identifier reference string

         DataValue getValue()   except + nogil  # wrap-doc:Returns the value of the term
         void setValue(DataValue value) except + nogil  # wrap-doc:Sets the value of the term

         void setUnit(Unit & unit) except + nogil  # wrap-doc:Sets the unit of the term
         Unit  getUnit() except + nogil  # wrap-doc:Returns the unit
         bool hasValue() except + nogil  # wrap-doc:Checks whether the term has a value
         bool hasUnit() except + nogil  # wrap-doc:Checks whether the term has a unit

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS::CVTerm":
    
    cdef cppclass Unit "OpenMS::CVTerm::Unit":
        Unit() except + nogil 
        Unit(Unit) except + nogil 
        String accession
        String name
        String cv_ref
        Unit(const String & p_accession, const String & p_name, const String & p_cv_ref) except + nogil 
        bool operator==(Unit & rhs) except + nogil 
        bool operator!=(Unit & rhs) except + nogil 
