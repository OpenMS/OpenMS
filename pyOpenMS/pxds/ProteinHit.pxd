from libcpp cimport bool
from Types cimport *
from DataValue cimport *
from Feature cimport *
from UniqueIdInterface cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/METADATA/ProteinHit.h>" namespace "OpenMS":

    cdef cppclass ProteinHit:


        ProteinHit() nogil except +
        ProteinHit(DoubleReal score, UInt rank, String accession, String sequence) nogil except +
        ProteinHit(ProteinHit) nogil except +

        Real getScore() nogil except +
        UInt getRank() nogil except +
        String getSequence() nogil except +
        String getAccession() nogil except +
        DoubleReal getCoverage() nogil except +

        void setScore(Real ) nogil except +
        void setRank(UInt) nogil except +
        void setSequence(String) nogil except +
        void setAccession(String) nogil except +
        void setCoverage(DoubleReal) nogil except +


        bool operator==(ProteinHit) nogil except +
        bool operator!=(ProteinHit) nogil except +
        bool isMetaEmpty() nogil except +
        void clearMetaInfo() nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys)
        void getKeys(libcpp_vector[unsigned int] & keys)
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +





