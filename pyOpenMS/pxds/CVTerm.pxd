from DataValue cimport *
from String cimport *
from libcpp cimport bool

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
