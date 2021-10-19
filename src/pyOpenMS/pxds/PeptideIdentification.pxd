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
        # wrap-doc:
        #   Represents the peptide hits for a spectrum
        #   -----
        #   This class is closely related to ProteinIdentification, which stores the protein hits
        #   and the general information about the identification run. More than one PeptideIdentification
        #   can belong to one ProteinIdentification. The general information about a
        #   PeptideIdentification has to be looked up in the corresponding ProteinIndentification, using
        #   the unique `identifier` that links the two.
        #   When loading PeptideHit instances from a File, the retention time and mass-to-charge ratio
        #   of the precursor spectrum can be accessed using getRT() and getMZ().
        #   This information can be used to map the peptide hits to an MSExperiment, a FeatureMap
        #   or a ConsensusMap using the IDMapper class

        PeptideIdentification() nogil except +
        PeptideIdentification(PeptideIdentification &) nogil except +
        bool operator==(PeptideIdentification) nogil except +
        bool operator!=(PeptideIdentification) nogil except +

        libcpp_vector[PeptideHit] getHits() nogil except + # wrap-doc:Returns the peptide hits as const
        void insertHit(PeptideHit) nogil except + # wrap-doc:Appends a peptide hit
        void setHits(libcpp_vector[PeptideHit]) nogil except + # wrap-doc:Sets the peptide hits

        double getSignificanceThreshold()   nogil except + # wrap-doc:Returns the peptide significance threshold value
        void setSignificanceThreshold(double value) nogil except + # wrap-doc:Setting of the peptide significance threshold value

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

        libcpp_vector[PeptideHit] getReferencingHits(libcpp_vector[PeptideHit], libcpp_set[String] &) nogil except + # wrap-doc:Returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)

