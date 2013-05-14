from libcpp cimport bool
from Types cimport *
from String cimport *
from DataValue cimport *

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS":

    cdef cppclass CVTerm:
         CVTerm()   nogil except +
         CVTerm(CVTerm)   nogil except +

         bool operator==(CVTerm)   nogil except +

         void setAccession(String accession) nogil except +
         String getAccession() nogil except +

         void setName(String name) nogil except +
         String getName() nogil except +

         void setCVIdentifierRef(String cv_id_ref) nogil except +
         String getCVIdentifierRef() nogil except +

         DataValue getValue()   nogil except +
         void setValue(DataValue value) nogil except +

         void setUnit(Unit & unit) nogil except +
         Unit  getUnit() nogil except +
         bool hasValue() nogil except +
         bool hasUnit() nogil except +

cdef extern from "<OpenMS/METADATA/CVTerm.h>" namespace "OpenMS::CVTerm":
    
    cdef cppclass Unit "OpenMS::CVTerm::Unit":
        Unit() nogil except +
        Unit(Unit) nogil except +
        String accession
        String name
        String cv_ref
        Unit(String & p_accession, String & p_name, String & p_cv_ref) nogil except +
        bool operator==(Unit & rhs) nogil except +
        bool operator!=(Unit & rhs) nogil except +

