from libcpp cimport bool
from Types cimport *
from DataValue cimport *
from Feature cimport *
from UniqueIdInterface cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ProteinHit.h>" namespace "OpenMS":

    cdef cppclass ProteinHit(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        ProteinHit() nogil except +
        ProteinHit(double score, UInt rank, String accession, String sequence) nogil except +
        ProteinHit(ProteinHit) nogil except +

        # const members
        ## double COVERAGE_UNKNOWN

        float getScore() nogil except +
        UInt getRank() nogil except +
        String getSequence() nogil except +
        String getAccession() nogil except +
        String getDescription() nogil except +
        double getCoverage() nogil except +

        void setScore(float ) nogil except +
        void setRank(UInt) nogil except +
        void setSequence(String) nogil except +
        void setAccession(String) nogil except +
        void setDescription(String description) nogil except +
        void setCoverage(double) nogil except +

        bool operator==(ProteinHit) nogil except +
        bool operator!=(ProteinHit) nogil except +



