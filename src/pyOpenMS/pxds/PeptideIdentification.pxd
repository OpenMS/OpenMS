from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoInterface cimport *
from PeptideHit cimport *

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

        bool       hasMZ() nogil except +
        double     getMZ() nogil except +
        void       setMZ(double) nogil except +

        bool       hasRT() nogil except +
        double     getRT() nogil except +
        void       setRT(double) nogil except +

        String     getBaseName() nogil except +
        void       setBaseName(String) nogil except +

        String     getExperimentLabel() nogil except +
        void       setExperimentLabel(String) nogil except +

        void       assignRanks() nogil except +
        void       sort() nogil except +
        void sortByRank() nogil except +
        bool       empty() nogil except +

        libcpp_vector[PeptideHit] getReferencingHits(libcpp_vector[PeptideHit], libcpp_set[String] &) nogil except +

