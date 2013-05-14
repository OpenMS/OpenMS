from libcpp cimport bool
from Types cimport *
from DataValue cimport *
from Feature cimport *
from ProteinIdentification cimport *
from AASequence cimport *

cdef extern from "<OpenMS/METADATA/PeptideHit.h>" namespace "OpenMS":

    cdef cppclass PeptideHit:


        PeptideHit() nogil except +

        PeptideHit(DoubleReal score,
                   UInt       rank,
                   Int        charge,
                   AASequence sequence
                   ) nogil except +

        PeptideHit(PeptideHit) nogil except + # wrap-ignore

        Real getScore() nogil except +
        UInt getRank() nogil except +
        AASequence getSequence() nogil except +
        libcpp_vector[String] getProteinAccessions() nogil except +
        void setProteinAccessions(libcpp_vector[String]) nogil except +


        void setScore(Real ) nogil except +
        void setRank(UInt) nogil except +
        void setSequence(AASequence) nogil except +
        void setCharge(Int) nogil except +

        void addProteinAccession(String) nogil except +
        void setAABefore(char) nogil except +
        char getAABefore() nogil except +
        void setAAAfter(char) nogil except +
        char getAAAfter() nogil except +


        bool operator==(PeptideHit) nogil except +
        bool operator!=(PeptideHit) nogil except +
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






