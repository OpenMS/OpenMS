from Types cimport *
from DataValue cimport *
from Feature cimport *
from AASequence cimport *
from PeptideEvidence cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/PeptideHit.h>" namespace "OpenMS":

    cdef cppclass PeptideHit(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface
        # wrap-doc:
        #  Represents a single peptide identification hit from a database search
        #  
        #  A PeptideHit stores information about a candidate peptide sequence that was
        #  matched to a spectrum. Each hit contains:
        #  - The peptide sequence (as AASequence)
        #  - A score from the search engine
        #  - The rank among all candidates
        #  - The charge state
        #  - Protein mappings (PeptideEvidence objects)
        #  
        #  Multiple PeptideHit objects are typically stored in a PeptideIdentification,
        #  sorted by score to show the most likely candidates first.
        #  
        #  Example usage:
        #    >>> hit = oms.PeptideHit()
        #    >>> hit.setSequence(oms.AASequence.fromString("PEPTIDER"))
        #    >>> hit.setScore(95.5)
        #    >>> hit.setRank(1)
        #    >>> hit.setCharge(2)
        #    >>> # Access information
        #    >>> print(f"Sequence: {hit.getSequence().toString()}")
        #    >>> print(f"Score: {hit.getScore()}, Rank: {hit.getRank()}")
        #    >>> print(f"Charge: {hit.getCharge()}")

        PeptideHit() except + nogil 

        PeptideHit(double score,
                   UInt       rank,
                   Int        charge,
                   AASequence sequence
                   ) except + nogil 

        PeptideHit(PeptideHit &) except + nogil 

        float getScore() except + nogil 
            # wrap-doc:
            #  Returns the score of this peptide-spectrum match (PSM)
            #  
            #  Returns:
            #    float: The search engine score
            #  
            #  Note:
            #    Interpretation depends on the score type (check isHigherScoreBetter)

        UInt getRank() except + nogil 
            # wrap-doc:
            #  Returns the rank of this hit among all candidates
            #  
            #  Returns:
            #    int: Rank (1 = best hit, 2 = second best, etc.)

        AASequence getSequence() except + nogil 
            # wrap-doc:
            #  Returns the peptide sequence
            #  
            #  Returns:
            #    AASequence: The peptide amino acid sequence with modifications

        Int getCharge() except + nogil 
            # wrap-doc:
            #  Returns the charge state of the peptide
            #  
            #  Returns:
            #    int: Charge state (e.g., 2 for doubly charged)

        libcpp_vector[PeptideEvidence] getPeptideEvidences() except + nogil 
            # wrap-doc:
            #  Returns protein mapping information for this peptide
            #  
            #  Returns:
            #    list of PeptideEvidence: List of proteins where this peptide was found
            #  
            #  Note:
            #    Each evidence contains protein accession, start/end positions, and if it's a decoy

        void setPeptideEvidences(libcpp_vector[PeptideEvidence]) except + nogil 
            # wrap-doc:
            #  Sets the protein mapping information
            #  
            #  Args:
            #    evidences (list of PeptideEvidence): Protein locations for this peptide

        void addPeptideEvidence(PeptideEvidence) except + nogil 
            # wrap-doc:
            #  Adds a single protein mapping
            #  
            #  Args:
            #    evidence (PeptideEvidence): Protein location information to add

        libcpp_set[String] extractProteinAccessionsSet() except + nogil 
            # wrap-doc:
            #  Extracts all unique protein accessions
            #  
            #  Returns:
            #    set of str: Set of unique protein accession strings
            #  
            #  Note:
            #    Empty accessions are excluded from the result

        bool isDecoy() except + nogil 
            # wrap-doc:
            #  Checks if this hit maps only to decoy proteins
            #  
            #  Returns:
            #    bool: True if all protein mappings are decoys, False otherwise
            #  
            #  Note:
            #    Returns False if no target/decoy information is available

        void setAnalysisResults(libcpp_vector[PeptideHit_AnalysisResult] aresult) except + nogil 
            # wrap-doc:
            #  Sets search engine sub-scores
            #  
            #  Args:
            #    aresult (list of PeptideHit_AnalysisResult): Sub-score information from search engine

        void addAnalysisResults(PeptideHit_AnalysisResult aresult) except + nogil 
            # wrap-doc:
            #  Adds a search engine sub-score
            #  
            #  Args:
            #    aresult (PeptideHit_AnalysisResult): Sub-score to add

        libcpp_vector[PeptideHit_AnalysisResult] getAnalysisResults() except + nogil 
            # wrap-doc:
            #  Returns all search engine sub-scores
            #  
            #  Returns:
            #    list of PeptideHit_AnalysisResult: Sub-score information

        void setPeakAnnotations(libcpp_vector[PeptideHit_PeakAnnotation]) except + nogil 
            # wrap-doc:
            #  Sets fragment ion annotations
            #  
            #  Args:
            #    annotations (list of PeptideHit_PeakAnnotation): Fragment peak annotations

        libcpp_vector[PeptideHit_PeakAnnotation] getPeakAnnotations() except + nogil 
            # wrap-doc:
            #  Returns fragment ion annotations
            #  
            #  Returns:
            #    list of PeptideHit_PeakAnnotation: Annotated fragment peaks

        void setScore(double) except + nogil 
            # wrap-doc:
            #  Sets the PSM score
            #  
            #  Args:
            #    score (float): The search engine score to set

        void setRank(UInt) except + nogil 
            # wrap-doc:
            #  Sets the rank of this hit
            #  
            #  Args:
            #    rank (int): Rank among all candidates (1 = best)

        void setSequence(AASequence) except + nogil 
            # wrap-doc:
            #  Sets the peptide sequence
            #  
            #  Args:
            #    sequence (AASequence): The peptide amino acid sequence

        void setCharge(Int) except + nogil 
            # wrap-doc:
            #  Sets the charge state
            #  
            #  Args:
            #    charge (int): Charge state of the peptide ion

        bool operator==(PeptideHit) except + nogil 
        bool operator!=(PeptideHit) except + nogil 

    cdef cppclass PeptideHit_AnalysisResult "OpenMS::PeptideHit::PepXMLAnalysisResult":

        PeptideHit_AnalysisResult() except + nogil  # compiler
        PeptideHit_AnalysisResult(PeptideHit_AnalysisResult &) except + nogil  # compiler

        String score_type
        bool higher_is_better
        double main_score
        libcpp_map[String, double] sub_scores


    cdef cppclass PeptideHit_PeakAnnotation "OpenMS::PeptideHit::PeakAnnotation":

        PeptideHit_PeakAnnotation() except + nogil  # compiler
        PeptideHit_PeakAnnotation(PeptideHit_PeakAnnotation &) except + nogil  # compiler

        String annotation
        int charge
        double mz
        double intensity

        void writePeakAnnotationsString_(String & annotation_string, libcpp_vector[ PeptideHit_PeakAnnotation ] annotations) except + nogil 

        bool operator<(PeptideHit_PeakAnnotation) except + nogil 
        bool operator==(PeptideHit_PeakAnnotation) except + nogil 

