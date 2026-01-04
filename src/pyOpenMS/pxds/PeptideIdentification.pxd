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
        #   MetaInfoInterface
        # wrap-hash:
        #  std
        # wrap-doc:
        #  Represents peptide identification results for a single spectrum or feature
        #  
        #  PeptideIdentification stores the results of peptide identification from database
        #  search engines (e.g., Mascot, X!Tandem, MSGF+). Each PeptideIdentification contains:
        #  
        #  - A list of peptide hits (candidate sequences) ranked by score
        #  - The precursor m/z and retention time
        #  - Score type and significance threshold
        #  - Link to the ProteinIdentification (via identifier)
        #  
        #  Multiple PeptideIdentifications can belong to one ProteinIdentification, which
        #  stores the search parameters and protein-level results.
        #  
        #  Example usage:
        #  
        #  .. code-block:: python
        #  
        #     pep_id = oms.PeptideIdentification()
        #     pep_id.setRT(1234.5)  # Set retention time
        #     pep_id.setMZ(445.678)  # Set precursor m/z
        #     pep_id.setScoreType("XTandem")
        #     # Add a peptide hit
        #     hit = oms.PeptideHit()
        #     hit.setScore(50.5)
        #     hit.setRank(1)
        #     hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
        #     hit.setCharge(2)
        #     pep_id.insertHit(hit)
        #     # Access hits
        #     for hit in pep_id.getHits():
        #         print(f"Sequence: {hit.getSequence().toString()}, Score: {hit.getScore()}")
        #  

        PeptideIdentification() except + nogil 
        PeptideIdentification(PeptideIdentification &) except + nogil 
        bool operator==(PeptideIdentification) except + nogil 
        bool operator!=(PeptideIdentification) except + nogil 

        libcpp_vector[PeptideHit] getHits() except + nogil 
            # wrap-doc:
            #  Returns all peptide hits (candidate sequences)
            #  
            #  :return: List of peptide candidates ranked by score
            #  
            #  Hits are typically sorted by score, with the best hit at index 0

        void insertHit(PeptideHit) except + nogil 
            # wrap-doc:
            #  Appends a peptide hit to the list
            #  
            #  :param hit: The peptide hit to add

        void setHits(libcpp_vector[PeptideHit]) except + nogil 
            # wrap-doc:
            #  Sets all peptide hits at once
            #  
            #  :param hits: List of peptide hits to store

        double getSignificanceThreshold()   except + nogil 
            # wrap-doc:
            #  Returns the significance threshold value
            #  
            #  :return: The threshold value (interpretation depends on score type)
            #  
            #  Hits with scores below/above this threshold (depending on score direction) may be considered insignificant

        void setSignificanceThreshold(double value) except + nogil 
            # wrap-doc:
            #  Sets the significance threshold value
            #  
            #  :param value: The threshold value to set

        String     getScoreType() except + nogil 
            # wrap-doc:
            #  Returns the type of score (e.g., "Mascot", "XTandem", "q-value")
            #  
            #  :return: Name of the score type

        void       setScoreType(String) except + nogil 
            # wrap-doc:
            #  Sets the score type
            #  
            #  :param score_type: Name of the score type (e.g., "Mascot", "XTandem")

        bool       isHigherScoreBetter() except + nogil 
            # wrap-doc:
            #  Returns whether higher scores are better
            #  
            #  :return: True if higher scores indicate better matches, False if lower is better

        void       setHigherScoreBetter(bool) except + nogil 
            # wrap-doc:
            #  Sets whether higher scores are better
            #  
            #  :param higher_better: True if higher scores are better, False otherwise

        String     getIdentifier() except + nogil 
            # wrap-doc:
            #  Returns the identifier linking to the parent ProteinIdentification
            #  
            #  :return: Unique identifier string
            #  
            #  Use this to find the corresponding ProteinIdentification with search parameters

        void       setIdentifier(String) except + nogil 
            # wrap-doc:
            #  Sets the identifier linking to a ProteinIdentification
            #  
            #  :param identifier: Unique identifier string

        bool       hasMZ() except + nogil 
            # wrap-doc:
            #  Checks if m/z value is set
            #  
            #  :return: True if m/z is available, False otherwise

        double     getMZ() except + nogil 
            # wrap-doc:
            #  Returns the precursor m/z value
            #  
            #  :return: Mass-to-charge ratio of the precursor ion

        void       setMZ(double) except + nogil 
            # wrap-doc:
            #  Sets the precursor m/z value
            #  
            #  :param mz: Mass-to-charge ratio of the precursor

        bool       hasRT() except + nogil 
            # wrap-doc:
            #  Checks if retention time is set
            #  
            #  :return: True if RT is available, False otherwise

        double     getRT() except + nogil 
            # wrap-doc:
            #  Returns the retention time of the precursor
            #  
            #  :return: Retention time in seconds

        void       setRT(double) except + nogil 
            # wrap-doc:
            #  Sets the retention time of the precursor
            #  
            #  :param rt: Retention time in seconds 

        String     getBaseName() except + nogil 
        void       setBaseName(String) except + nogil 

        String     getExperimentLabel() except + nogil 
        void       setExperimentLabel(String) except + nogil 

        void       sort() except + nogil 
        bool       empty() except + nogil 

        libcpp_vector[PeptideHit] getReferencingHits(libcpp_vector[PeptideHit], libcpp_set[String] &) except + nogil  # wrap-doc:Returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)

 

