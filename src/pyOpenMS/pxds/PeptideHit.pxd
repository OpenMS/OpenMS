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

        PeptideHit() nogil except +

        PeptideHit(double score,
                   UInt       rank,
                   Int        charge,
                   AASequence sequence
                   ) nogil except +

        PeptideHit(PeptideHit) nogil except + # wrap-ignore

        float getScore() nogil except +
        UInt getRank() nogil except +
        AASequence getSequence() nogil except +
        Int getCharge() nogil except +
        libcpp_vector[PeptideEvidence] getPeptideEvidences() nogil except +
        void setPeptideEvidences(libcpp_vector[PeptideEvidence]) nogil except +
        void addPeptideEvidence(PeptideEvidence) nogil except +
        libcpp_set[String] extractProteinAccessionsSet() nogil except +

        void setAnalysisResults(libcpp_vector[PeptideHit_AnalysisResult] aresult) nogil except +
        void addAnalysisResults(PeptideHit_AnalysisResult aresult) nogil except +
        libcpp_vector[PeptideHit_AnalysisResult] getAnalysisResults() nogil except +

        void setPeakAnnotations(libcpp_vector[PeptideHit_PeakAnnotation]) nogil except +
        libcpp_vector[PeptideHit_PeakAnnotation] getPeakAnnotations() nogil except +

        void setScore(double) nogil except +
        void setRank(UInt) nogil except +
        void setSequence(AASequence) nogil except +
        void setCharge(Int) nogil except +

        bool operator==(PeptideHit) nogil except +
        bool operator!=(PeptideHit) nogil except +

    cdef cppclass PeptideHit_AnalysisResult "OpenMS::PeptideHit::PepXMLAnalysisResult":

        PeptideHit_AnalysisResult() nogil except +
        PeptideHit_AnalysisResult(PeptideHit_AnalysisResult) nogil except + #wrap-ignore

        String score_type
        bool higher_is_better
        double main_score
        libcpp_map[String, double] sub_scores


    cdef cppclass PeptideHit_PeakAnnotation "OpenMS::PeptideHit::PeakAnnotation":

        PeptideHit_PeakAnnotation() nogil except +
        PeptideHit_PeakAnnotation(PeptideHit_PeakAnnotation) nogil except + # wrap-ignore

        String annotation
        int charge
        double mz
        double intensity

        void writePeakAnnotationsString_(String & annotation_string, libcpp_vector[ PeptideHit_PeakAnnotation ] annotations) nogil except +

        bool operator<(PeptideHit_PeakAnnotation) nogil except +
        bool operator==(PeptideHit_PeakAnnotation) nogil except +

