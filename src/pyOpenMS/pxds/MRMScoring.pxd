from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from ISpectrumAccess cimport *
# from ITransition cimport *
# from TransitionExperiment cimport *
# from StatsHelpers cimport *
# from Scoring cimport *
from LightTargetedExperiment cimport *
from Matrix cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>" namespace "OpenSwath":
    
    cdef cppclass MRMScoring:
        MRMScoring() nogil except + # compiler
        MRMScoring(MRMScoring &) nogil except + # compiler

        # TODO create class for XCorrMatrix
        # XCorrMatrixType  getXCorrMatrix() nogil except +
        # NAMESPACE # # POINTER # void initializeXCorrMatrix(OpenSwath::IMRMFeature * mrmfeature, OpenSwath::ITransitionGroup * transition_group, bool normalize) nogil except +
        double calcXcorrCoelutionScore() nogil except + # wrap-doc:Calculate the cross-correlation coelution score. The score is a distance where zero indicates perfect coelution
        double calcXcorrCoelutionWeightedScore(libcpp_vector[ double ] & normalized_library_intensity) nogil except +
            # wrap-doc:
                #   Calculate the weighted cross-correlation coelution score
                #   -----
                #   The score is a distance where zero indicates perfect coelution. The
                #   score is weighted by the transition intensities, non-perfect coelution
                #   in low-intensity transitions should thus become less important

        libcpp_vector[ double ] calcSeparateXcorrContrastCoelutionScore() nogil except + # wrap-doc:Calculate the separate cross-correlation contrast score
        double calcXcorrPrecursorContrastCoelutionScore() nogil except +
            # wrap-doc:
                #   Calculate the precursor cross-correlation contrast score against the transitions
                #   -----
                #   The score is a distance where zero indicates perfect coelution

        double calcXcorrShapeScore() nogil except +
            # wrap-doc:
                #   Calculate the cross-correlation shape score
                #   -----
                #   The score is a correlation measure where 1 indicates perfect correlation
                #   and 0 means no correlation.

        double calcXcorrShapeWeightedScore(libcpp_vector[ double ] & normalized_library_intensity) nogil except +
            # wrap-doc:
                #   Calculate the weighted cross-correlation shape score
                #   -----
                #   The score is a correlation measure where 1 indicates perfect correlation
                #   and 0 means no correlation. The score is weighted by the transition
                #   intensities, non-perfect coelution in low-intensity transitions should
                #   thus become less important

        libcpp_vector[ double ] calcSeparateXcorrContrastShapeScore() nogil except + # wrap-doc:Calculate the separate cross-correlation contrast shape score
        double calcXcorrPrecursorContrastShapeScore() nogil except + # wrap-doc:Calculate the precursor cross-correlation shape score against the transitions

        # NAMESPACE # # POINTER # void calcLibraryScore(OpenSwath::IMRMFeature * mrmfeature, libcpp_vector[ TransitionType ] & transitions, double & correlation, double & rmsd, double & manhattan, double & dotprod) nogil except +
        double calcRTScore(LightCompound & peptide, double normalized_experimental_rt) nogil except +
        # NAMESPACE # # POINTER # double calcSNScore(OpenSwath::IMRMFeature * mrmfeature, libcpp_vector[ OpenSwath::ISignalToNoisePtr ] & signal_noise_estimators) nogil except +

        double calcMIScore() nogil except +
        double calcMIWeightedScore(const libcpp_vector[ double ] & normalized_library_intensity) nogil except +
        
        double calcMIPrecursorScore() nogil except + 
        double calcMIPrecursorContrastScore() nogil except + 
        double calcMIPrecursorCombinedScore() nogil except + 
        libcpp_vector[ double ] calcSeparateMIContrastScore() nogil except +              
        
        Matrix[ double ] getMIMatrix() nogil except +
