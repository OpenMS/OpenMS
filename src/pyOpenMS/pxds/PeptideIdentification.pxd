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
        # wrap-doc:
        #  Represents the peptide hits for a spectrum
        #  
        #  This class is closely related to ProteinIdentification, which stores the protein hits
        #  and the general information about the identification run. More than one PeptideIdentification
        #  can belong to one ProteinIdentification. The general information about a
        #  PeptideIdentification has to be looked up in the corresponding ProteinIndentification, using
        #  the unique `identifier` that links the two.
        #  When loading PeptideHit instances from a File, the retention time and mass-to-charge ratio
        #  of the precursor spectrum can be accessed using getRT() and getMZ().
        #  This information can be used to map the peptide hits to an MSExperiment, a FeatureMap
        #  or a ConsensusMap using the IDMapper class

        PeptideIdentification() except + nogil 
        PeptideIdentification(PeptideIdentification &) except + nogil 
        bool operator==(PeptideIdentification) except + nogil 
        bool operator!=(PeptideIdentification) except + nogil 

        libcpp_vector[PeptideHit] getHits() except + nogil  # wrap-doc:Returns the peptide hits as const
        void insertHit(PeptideHit) except + nogil  # wrap-doc:Appends a peptide hit
        void setHits(libcpp_vector[PeptideHit]) except + nogil  # wrap-doc:Sets the peptide hits

        double getSignificanceThreshold()   except + nogil  # wrap-doc:Returns the peptide significance threshold value
        void setSignificanceThreshold(double value) except + nogil  # wrap-doc:Setting of the peptide significance threshold value

        String     getScoreType() except + nogil 
        void       setScoreType(String) except + nogil 
        bool       isHigherScoreBetter() except + nogil 
        void       setHigherScoreBetter(bool) except + nogil 
        String     getIdentifier() except + nogil 
        void       setIdentifier(String) except + nogil 

        bool       hasMZ() except + nogil 
        double     getMZ() except + nogil 
        void       setMZ(double) except + nogil 

        bool       hasRT() except + nogil 
        double     getRT() except + nogil 
        void       setRT(double) except + nogil 

        String     getBaseName() except + nogil 
        void       setBaseName(String) except + nogil 

        String     getExperimentLabel() except + nogil 
        void       setExperimentLabel(String) except + nogil 

        void       assignRanks() except + nogil 
        void       sort() except + nogil 
        void sortByRank() except + nogil 
        bool       empty() except + nogil 

        libcpp_vector[PeptideHit] getReferencingHits(libcpp_vector[PeptideHit], libcpp_set[String] &) except + nogil  # wrap-doc:Returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)

