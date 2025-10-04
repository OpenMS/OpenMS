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
        #  Represents a single protein identification hit from a database search
        #  
        #  A ProteinHit stores information about a protein that was identified based on
        #  peptide evidence. Each hit contains:
        #  - Protein accession (database identifier)
        #  - Score from protein inference
        #  - Rank among protein candidates
        #  - Protein sequence (optional)
        #  - Sequence coverage percentage
        #  
        #  Multiple ProteinHit objects are stored in a ProteinIdentification, typically
        #  sorted by score to show the most confident identifications first.
        #  
        #  Example usage:
        #    >>> protein_hit = oms.ProteinHit()
        #    >>> protein_hit.setAccession("P12345")
        #    >>> protein_hit.setScore(150.5)
        #    >>> protein_hit.setRank(1)
        #    >>> protein_hit.setCoverage(45.2)  # 45.2% coverage
        #    >>> protein_hit.setDescription("Example protein")
        #    >>> # Access information
        #    >>> print(f"Accession: {protein_hit.getAccession()}")
        #    >>> print(f"Score: {protein_hit.getScore()}, Coverage: {protein_hit.getCoverage()}%")

        ProteinHit() except + nogil
        ProteinHit(double score, UInt rank, String accession, String sequence) except + nogil 
        ProteinHit(ProteinHit &) except + nogil 

        # const members
        ## double COVERAGE_UNKNOWN

        float getScore() except + nogil 
            # wrap-doc:
            #  Returns the protein inference score
            #  
            #  Returns:
            #    float: Score from protein inference algorithm

        UInt getRank() except + nogil 
            # wrap-doc:
            #  Returns the rank of this protein hit
            #  
            #  Returns:
            #    int: Rank (1 = best hit, 2 = second best, etc.)

        String getSequence() except + nogil 
            # wrap-doc:
            #  Returns the protein sequence
            #  
            #  Returns:
            #    str: Full amino acid sequence of the protein (if available)

        String getAccession() except + nogil 
            # wrap-doc:
            #  Returns the protein accession
            #  
            #  Returns:
            #    str: Database accession/identifier (e.g., "P12345" for UniProt)

        String getDescription() except + nogil 
            # wrap-doc:
            #  Returns the protein description
            #  
            #  Returns:
            #    str: Human-readable protein name/description from database

        double getCoverage() except + nogil 
            # wrap-doc:
            #  Returns the sequence coverage percentage
            #  
            #  Returns:
            #    float: Percentage of protein sequence covered by identified peptides
            #  
            #  Note:
            #    Value is in range 0-100 (e.g., 45.2 means 45.2% coverage)

        void setScore(float ) except + nogil 
            # wrap-doc:
            #  Sets the protein inference score
            #  
            #  Args:
            #    score (float): Score to set

        void setRank(UInt) except + nogil 
            # wrap-doc:
            #  Sets the rank
            #  
            #  Args:
            #    rank (int): Rank among all protein candidates

        void setSequence(String) except + nogil 
            # wrap-doc:
            #  Sets the protein sequence
            #  
            #  Args:
            #    sequence (str): Full amino acid sequence

        void setAccession(String) except + nogil 
            # wrap-doc:
            #  Sets the protein accession
            #  
            #  Args:
            #    accession (str): Database accession/identifier

        void setDescription(String description) except + nogil 
            # wrap-doc:
            #  Sets the protein description
            #  
            #  Args:
            #    description (str): Human-readable protein name/description

        void setCoverage(double) except + nogil 
            # wrap-doc:
            #  Sets the sequence coverage percentage
            #  
            #  Args:
            #    coverage (float): Percentage (0-100) of sequence covered by peptides

        bool isDecoy() except + nogil 
            # wrap-doc:
            #  Checks if this is a decoy protein hit
            #  
            #  Returns:
            #    bool: True if this is a decoy hit from target-decoy search

        bool operator==(ProteinHit) except + nogil
        bool operator!=(ProteinHit) except + nogil
