from Types cimport *
from DataValue cimport *
from Feature cimport *
from AASequence cimport *
from PeptideEvidence cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/PeptideHit.h>" namespace "OpenMS":

    cdef cppclass PeptideHit(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface
        # wrap-doc:
        #   Representation of a peptide hit
        #   -----
        #   It contains the fields score, score_type, rank, and sequence

        PeptideHit() nogil except +

        PeptideHit(double score,
                   UInt       rank,
                   Int        charge,
                   AASequence sequence
                   ) nogil except +

        PeptideHit(PeptideHit &) nogil except +

        float getScore() nogil except + Returns the PSM score
        UInt getRank() nogil except + Returns the PSM rank
        AASequence getSequence() nogil except + Returns the peptide sequence without trailing or following spaces
        Int getCharge() nogil except + Returns the charge of the peptide
        libcpp_vector[PeptideEvidence] getPeptideEvidences() nogil except + Returns information on peptides (potentially) identified by this PSM
        void setPeptideEvidences(libcpp_vector[PeptideEvidence]) nogil except + Set information on peptides (potentially) identified by this PSM
        void addPeptideEvidence(PeptideEvidence) nogil except + Adds information on a peptide that is (potentially) identified by this PSM
        libcpp_set[String] extractProteinAccessionsSet() nogil except + Extracts the set of non-empty protein accessions from peptide evidences

        void setAnalysisResults(libcpp_vector[PeptideHit_AnalysisResult] aresult) nogil except + Set information on (search engine) sub scores associated with this PSM
        void addAnalysisResults(PeptideHit_AnalysisResult aresult) nogil except + Add information on (search engine) sub scores associated with this PSM
        libcpp_vector[PeptideHit_AnalysisResult] getAnalysisResults() nogil except + Returns information on (search engine) sub scores associated with this PSM

        void setPeakAnnotations(libcpp_vector[PeptideHit_PeakAnnotation]) nogil except + Sets the fragment annotations
        libcpp_vector[PeptideHit_PeakAnnotation] getPeakAnnotations() nogil except +

        void setScore(double) nogil except + Sets the PSM score
        void setRank(UInt) nogil except + Sets the PSM rank
        void setSequence(AASequence) nogil except + Sets the peptide sequence
        void setCharge(Int) nogil except + Sets the charge of the peptide

        bool operator==(PeptideHit) nogil except +
        bool operator!=(PeptideHit) nogil except +

    cdef cppclass PeptideHit_AnalysisResult "OpenMS::PeptideHit::PepXMLAnalysisResult":

        PeptideHit_AnalysisResult() nogil except + # compiler
        PeptideHit_AnalysisResult(PeptideHit_AnalysisResult &) nogil except + # compiler

        String score_type
        bool higher_is_better
        double main_score
        libcpp_map[String, double] sub_scores


    cdef cppclass PeptideHit_PeakAnnotation "OpenMS::PeptideHit::PeakAnnotation":

        PeptideHit_PeakAnnotation() nogil except + # compiler
        PeptideHit_PeakAnnotation(PeptideHit_PeakAnnotation &) nogil except + # compiler

        String annotation
        int charge
        double mz
        double intensity

        void writePeakAnnotationsString_(String & annotation_string, libcpp_vector[ PeptideHit_PeakAnnotation ] annotations) nogil except +

        bool operator<(PeptideHit_PeakAnnotation) nogil except +
        bool operator==(PeptideHit_PeakAnnotation) nogil except +

