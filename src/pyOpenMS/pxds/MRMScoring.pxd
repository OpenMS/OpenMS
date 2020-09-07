from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from ISpectrumAccess cimport *
# from ITransition cimport *
# from TransitionExperiment cimport *
# from StatsHelpers cimport *
# from Scoring cimport *
from LightTargetedExperiment cimport *

cdef extern from "<OpenMS/OPENSWATHALGO/ALGO/MRMScoring.h>" namespace "OpenSwath":
    
    cdef cppclass MRMScoring "OpenSwath::MRMScoring":
        MRMScoring(MRMScoring) nogil except + #wrap-ignore

        # TODO create class for XCorrMatrix
        # XCorrMatrixType  getXCorrMatrix() nogil except +
        # NAMESPACE # # POINTER # void initializeXCorrMatrix(OpenSwath::IMRMFeature * mrmfeature, OpenSwath::ITransitionGroup * transition_group, bool normalize) nogil except +
        double calcXcorrCoelutionScore() nogil except +
        double calcXcorrCoelutionWeightedScore(libcpp_vector[ double ] & normalized_library_intensity) nogil except +
        libcpp_vector[ double ] calcSeparateXcorrContrastCoelutionScore() nogil except +
        double calcXcorrPrecursorContrastCoelutionScore() nogil except +

        double calcXcorrShapeScore() nogil except +
        double calcXcorrShapeWeightedScore(libcpp_vector[ double ] & normalized_library_intensity) nogil except +
        libcpp_vector[ double ] calcSeparateXcorrContrastShapeScore() nogil except +
        double calcXcorrPrecursorContrastShapeScore() nogil except +

        # NAMESPACE # # POINTER # void calcLibraryScore(OpenSwath::IMRMFeature * mrmfeature, libcpp_vector[ TransitionType ] & transitions, double & correlation, double & rmsd, double & manhattan, double & dotprod) nogil except +
        double calcRTScore(LightCompound & peptide, double normalized_experimental_rt) nogil except +
        # NAMESPACE # # POINTER # double calcSNScore(OpenSwath::IMRMFeature * mrmfeature, libcpp_vector[ OpenSwath::ISignalToNoisePtr ] & signal_noise_estimators) nogil except +

        double calcMIScore() nogil except +
        double calcMIWeightedScore(const libcpp_vector[ double ] & normalized_library_intensity) nogil except +
        
        double calcMIPrecursorScore() nogil except + 
        double calcMIPrecursorContrastScore() nogil except + 
        double calcMIPrecursorCombinedScore() nogil except + 
        libcpp_vector[ double ] calcSeparateMIContrastScore() nogil except +              
        
        libcpp_vector[ libcpp_vector[ double ] ]  getMIMatrix() nogil except +
