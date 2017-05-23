from Types cimport *
from DataValue cimport *
from Feature cimport *
from ProteinIdentification cimport *
from AASequence cimport *
from PeptideEvidence cimport *

cdef extern from "<OpenMS/METADATA/PeptideHit.h>" namespace "OpenMS":

    cdef cppclass PeptideHit:


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

        void setFragmentAnnotations(libcpp_vector[PeptideHit_FragmentAnnotation]) nogil except +
        libcpp_vector[PeptideHit_FragmentAnnotation] getFragmentAnnotations() nogil except +

        void setScore(double) nogil except +
        void setRank(UInt) nogil except +
        void setSequence(AASequence) nogil except +
        void setCharge(Int) nogil except +

        bool operator==(PeptideHit) nogil except +
        bool operator!=(PeptideHit) nogil except +
        bool isMetaEmpty() nogil except +
        void clearMetaInfo() nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +


    cdef cppclass PeptideHit_AnalysisResult "OpenMS::PeptideHit::PepXMLAnalysisResult":

        PeptideHit_AnalysisResult() nogil except +
        PeptideHit_AnalysisResult(PeptideHit_AnalysisResult) nogil except + #wrap-ignore

        String score_type
        bool higher_is_better
        double main_score
        libcpp_map[String, double] sub_scores


    cdef cppclass PeptideHit_FragmentAnnotation "OpenMS::PeptideHit::FragmentAnnotation":

        PeptideHit_FragmentAnnotation() nogil except +
        PeptideHit_FragmentAnnotation(PeptideHit_FragmentAnnotation) nogil except + # wrap-ignore

        String annotation
        int charge
        double mz
        double intensity

        bool operator<(PeptideHit_FragmentAnnotation) nogil except +
        bool operator==(PeptideHit_FragmentAnnotation) nogil except +

