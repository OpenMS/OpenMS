from libcpp cimport bool
from Types cimport *
from String cimport *
from DataValue cimport *

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS":

    cdef cppclass CVTerm:
         CVTerm() nogil except +
         CVTerm(CVTerm &) nogil except +

         bool operator==(CVTerm) nogil except +

         void setAccession(String accession) nogil except + # wrap-doc:Sets the accession string of the term
         String getAccession() nogil except + # wrap-doc:Returns the accession string of the term

         void setName(String name) nogil except + # wrap-doc:Sets the name of the term
         String getName() nogil except + # wrap-doc:Returns the name of the term

         void setCVIdentifierRef(String cv_id_ref) nogil except + # wrap-doc:Sets the CV identifier reference string, e.g. UO for unit obo
         String getCVIdentifierRef() nogil except + # wrap-doc:Returns the CV identifier reference string

         DataValue getValue()   nogil except + # wrap-doc:Returns the value of the term
         void setValue(DataValue value) nogil except + # wrap-doc:Sets the value of the term

         void setUnit(Unit & unit) nogil except + # wrap-doc:Sets the unit of the term
         Unit  getUnit() nogil except + # wrap-doc:Returns the unit
         bool hasValue() nogil except + # wrap-doc:Checks whether the term has a value
         bool hasUnit() nogil except + # wrap-doc:Checks whether the term has a unit

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS::CVTerm":
    
    cdef cppclass Unit "OpenMS::CVTerm::Unit":
        Unit() nogil except +
        Unit(Unit) nogil except +
        String accession
        String name
        String cv_ref
        Unit(const String & p_accession, const String & p_name, const String & p_cv_ref) nogil except +
        bool operator==(Unit & rhs) nogil except +
        bool operator!=(Unit & rhs) nogil except +
