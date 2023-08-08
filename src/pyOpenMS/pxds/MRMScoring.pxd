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
        MRMScoring() except + nogil  # compiler
        MRMScoring(MRMScoring &) except + nogil  # compiler

        # TODO create class for XCorrMatrix
        # XCorrMatrixType  getXCorrMatrix() except + nogil 
        # NAMESPACE # # POINTER # void initializeXCorrMatrix(OpenSwath::IMRMFeature * mrmfeature, OpenSwath::ITransitionGroup * transition_group, bool normalize) except + nogil 
        double calcXcorrCoelutionScore() except + nogil  # wrap-doc:Calculate the cross-correlation coelution score. The score is a distance where zero indicates perfect coelution
        double calcXcorrCoelutionWeightedScore(libcpp_vector[ double ] & normalized_library_intensity) except + nogil 
            # wrap-doc:
                #  Calculate the weighted cross-correlation coelution score
                #  
                #  The score is a distance where zero indicates perfect coelution. The
                #  score is weighted by the transition intensities, non-perfect coelution
                #  in low-intensity transitions should thus become less important

        libcpp_vector[ double ] calcSeparateXcorrContrastCoelutionScore() except + nogil  # wrap-doc:Calculate the separate cross-correlation contrast score
        double calcXcorrPrecursorContrastCoelutionScore() except + nogil 
            # wrap-doc:
                #  Calculate the precursor cross-correlation contrast score against the transitions
                #  
                #  The score is a distance where zero indicates perfect coelution

        double calcXcorrShapeScore() except + nogil 
            # wrap-doc:
                #  Calculate the cross-correlation shape score
                #  
                #  The score is a correlation measure where 1 indicates perfect correlation
                #  and 0 means no correlation.

        double calcXcorrShapeWeightedScore(libcpp_vector[ double ] & normalized_library_intensity) except + nogil 
            # wrap-doc:
                #  Calculate the weighted cross-correlation shape score
                #  
                #  The score is a correlation measure where 1 indicates perfect correlation
                #  and 0 means no correlation. The score is weighted by the transition
                #  intensities, non-perfect coelution in low-intensity transitions should
                #  thus become less important

        libcpp_vector[ double ] calcSeparateXcorrContrastShapeScore() except + nogil  # wrap-doc:Calculate the separate cross-correlation contrast shape score
        double calcXcorrPrecursorContrastShapeScore() except + nogil  # wrap-doc:Calculate the precursor cross-correlation shape score against the transitions

        # NAMESPACE # # POINTER # void calcLibraryScore(OpenSwath::IMRMFeature * mrmfeature, libcpp_vector[ TransitionType ] & transitions, double & correlation, double & rmsd, double & manhattan, double & dotprod) except + nogil 
        double calcRTScore(LightCompound & peptide, double normalized_experimental_rt) except + nogil 
        # NAMESPACE # # POINTER # double calcSNScore(OpenSwath::IMRMFeature * mrmfeature, libcpp_vector[ OpenSwath::ISignalToNoisePtr ] & signal_noise_estimators) except + nogil 

        double calcMIScore() except + nogil 
        double calcMIWeightedScore(const libcpp_vector[ double ] & normalized_library_intensity) except + nogil 
        
        double calcMIPrecursorScore() except + nogil  
        double calcMIPrecursorContrastScore() except + nogil  
        double calcMIPrecursorCombinedScore() except + nogil  
        libcpp_vector[ double ] calcSeparateMIContrastScore() except + nogil               
        
        Matrix[ double ] getMIMatrix() except + nogil 
