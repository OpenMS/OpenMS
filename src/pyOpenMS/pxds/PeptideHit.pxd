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
        #  Representation of a peptide hit
        #  
        #  It contains the fields score, score_type, rank, and sequence

        PeptideHit() except + nogil 

        PeptideHit(double score,
                   UInt       rank,
                   Int        charge,
                   AASequence sequence
                   ) except + nogil 

        PeptideHit(PeptideHit &) except + nogil 

        float getScore() except + nogil  # wrap-doc:Returns the PSM score
        UInt getRank() except + nogil  # wrap-doc:Returns the PSM rank
        AASequence getSequence() except + nogil  # wrap-doc:Returns the peptide sequence without trailing or following spaces
        Int getCharge() except + nogil  # wrap-doc:Returns the charge of the peptide
        libcpp_vector[PeptideEvidence] getPeptideEvidences() except + nogil  # wrap-doc:Returns information on peptides (potentially) identified by this PSM
        void setPeptideEvidences(libcpp_vector[PeptideEvidence]) except + nogil  # wrap-doc:Sets information on peptides (potentially) identified by this PSM
        void addPeptideEvidence(PeptideEvidence) except + nogil  # wrap-doc:Adds information on a peptide that is (potentially) identified by this PSM
        libcpp_set[String] extractProteinAccessionsSet() except + nogil  # wrap-doc:Extracts the set of non-empty protein accessions from peptide evidences

        void setAnalysisResults(libcpp_vector[PeptideHit_AnalysisResult] aresult) except + nogil  # wrap-doc:Sets information on (search engine) sub scores associated with this PSM
        void addAnalysisResults(PeptideHit_AnalysisResult aresult) except + nogil  # wrap-doc:Add information on (search engine) sub scores associated with this PSM
        libcpp_vector[PeptideHit_AnalysisResult] getAnalysisResults() except + nogil  # wrap-doc:Returns information on (search engine) sub scores associated with this PSM

        void setPeakAnnotations(libcpp_vector[PeptideHit_PeakAnnotation]) except + nogil  # wrap-doc:Sets the fragment annotations
        libcpp_vector[PeptideHit_PeakAnnotation] getPeakAnnotations() except + nogil  # wrap-doc:Returns the fragment annotations

        void setScore(double) except + nogil  # wrap-doc:Sets the PSM score
        void setRank(UInt) except + nogil  # wrap-doc:Sets the PSM rank
        void setSequence(AASequence) except + nogil  # wrap-doc:Sets the peptide sequence
        void setCharge(Int) except + nogil  # wrap-doc:Sets the charge of the peptide

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

