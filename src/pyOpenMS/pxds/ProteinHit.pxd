from libcpp cimport bool
from Types cimport *
from DataValue cimport *
from Feature cimport *
from UniqueIdInterface cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/ProteinHit.h>" namespace "OpenMS":

    cdef cppclass ProteinHit(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface
        # wrap-doc:
        #  Representation of a protein hit
        #  
        #  It contains the fields score, score_type, rank, accession,
        #  sequence and coverage

        ProteinHit() except + nogil 
        ProteinHit(double score, UInt rank, String accession, String sequence) except + nogil 
        ProteinHit(ProteinHit &) except + nogil 

        # const members
        ## double COVERAGE_UNKNOWN

        float getScore() except + nogil  # wrap-doc:Returns the score of the protein hit
        UInt getRank() except + nogil  # wrap-doc:Returns the rank of the protein hit
        String getSequence() except + nogil  # wrap-doc:Returns the protein sequence
        String getAccession() except + nogil  # wrap-doc:Returns the accession of the protein
        String getDescription() except + nogil  # wrap-doc:Returns the description of the protein
        double getCoverage() except + nogil  # wrap-doc:Returns the coverage (in percent) of the protein hit based upon matched peptides

        void setScore(float ) except + nogil  # wrap-doc:Sets the score of the protein hit
        void setRank(UInt) except + nogil  # wrap-doc:Sets the rank
        void setSequence(String) except + nogil  # wrap-doc:Sets the protein sequence
        void setAccession(String) except + nogil  # wrap-doc:Sets the accession of the protein
        void setDescription(String description) except + nogil  # wrap-doc:Sets the description of the protein
        void setCoverage(double) except + nogil  # wrap-doc:Sets the coverage (in percent) of the protein hit based upon matched peptides

        bool operator==(ProteinHit) except + nogil 
        bool operator!=(ProteinHit) except + nogil 
