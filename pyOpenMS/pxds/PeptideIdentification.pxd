from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoInterface cimport *
from PeptideHit cimport *
from ProteinHit cimport *

cdef extern from "<OpenMS/METADATA/PeptideIdentification.h>" namespace "OpenMS":

    cdef cppclass PeptideIdentification(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        PeptideIdentification() nogil except +
        PeptideIdentification(PeptideIdentification) nogil except +
        bool operator==(PeptideIdentification) nogil except +
        bool operator!=(PeptideIdentification) nogil except +

        libcpp_vector[PeptideHit] getHits() nogil except +
        void insertHit(PeptideHit) nogil except +
        void setHits(libcpp_vector[PeptideHit]) nogil except +

        double getSignificanceThreshold()   nogil except +
        # setting of the peptide significance threshold value
        void setSignificanceThreshold(double value) nogil except +

        String     getScoreType() nogil except +
        void       setScoreType(String) nogil except +
        bool       isHigherScoreBetter() nogil except +
        void       setHigherScoreBetter(bool) nogil except +
        String     getIdentifier() nogil except +
        void       setIdentifier(String) nogil except +

        String     getBaseName() nogil except +
        void       setBaseName(String) nogil except +

        void       assignRanks() nogil except +
        void       sort() nogil except +
        bool       empty() nogil except +

        void       getReferencingHits(String, libcpp_vector[PeptideHit] &) nogil except +
        void       getReferencingHits(libcpp_vector[String], libcpp_vector[PeptideHit] &) nogil except +
        void       getReferencingHits(libcpp_vector[ProteinHit], libcpp_vector[PeptideHit] &) nogil except +

        void       getNonReferencingHits(String, libcpp_vector[PeptideHit] &) nogil except +
        void       getNonReferencingHits(libcpp_vector[String], libcpp_vector[PeptideHit] &) nogil except +
        void       getNonReferencingHits(libcpp_vector[ProteinHit], libcpp_vector[PeptideHit] &) nogil except +


        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
