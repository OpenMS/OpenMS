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
        # wrap-doc:
        #   Representation of a protein hit
        #   -----
        #   It contains the fields score, score_type, rank, accession,
        #   sequence and coverage

        ProteinHit() nogil except +
        ProteinHit(double score, UInt rank, String accession, String sequence) nogil except +
        ProteinHit(ProteinHit &) nogil except +

        # const members
        ## double COVERAGE_UNKNOWN

        float getScore() nogil except + # wrap-doc:Returns the score of the protein hit
        UInt getRank() nogil except + # wrap-doc:Returns the rank of the protein hit
        String getSequence() nogil except + # wrap-doc:Returns the protein sequence
        String getAccession() nogil except + # wrap-doc:Returns the accession of the protein
        String getDescription() nogil except + # wrap-doc:Returns the description of the protein
        double getCoverage() nogil except + # wrap-doc:Returns the coverage (in percent) of the protein hit based upon matched peptides

        void setScore(float ) nogil except + # wrap-doc:Sets the score of the protein hit
        void setRank(UInt) nogil except + # wrap-doc:Sets the rank
        void setSequence(String) nogil except + # wrap-doc:Sets the protein sequence
        void setAccession(String) nogil except + # wrap-doc:Sets the accession of the protein
        void setDescription(String description) nogil except + # wrap-doc:Sets the description of the protein
        void setCoverage(double) nogil except + # wrap-doc:Sets the coverage (in percent) of the protein hit based upon matched peptides

        bool operator==(ProteinHit) nogil except +
        bool operator!=(ProteinHit) nogil except +
