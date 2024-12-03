from smart_ptr cimport shared_ptr
from Types cimport *
from OpenSwathDataStructures cimport *
from SpectrumAccessOpenMS cimport *
from LightTargetedExperiment cimport LightTransition

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>" namespace "OpenMS":

    cdef cppclass OpenSwathScoring:

        OpenSwathScoring() except + nogil 
        OpenSwathScoring(OpenSwathScoring &) except + nogil  # compiler

        void initialize(double rt_normalization_factor,
                        int add_up_spectra,
                        double spacing_for_spectra_resampling,
                        double drift_extra,
                        OpenSwath_Scores_Usage su,
                        libcpp_string spectrum_addition_method,
                        bool use_ms1_ion_mobility) except + nogil
            # wrap-doc:
                #  Initialize the scoring object\n
                #  Sets the parameters for the scoring
                #  
                #  
                #  :param rt_normalization_factor: Specifies the range of the normalized retention time space
                #  :param add_up_spectra: How many spectra to add up (default 1)
                #  :param spacing_for_spectra_resampling: Spacing factor for spectra addition
                #  :param su: Which scores to actually compute
                #  :param spectrum_addition_method: Method to use for spectrum addition (valid: "simple", "resample")

        # void calculateChromatographicScores(
        #      OpenSwath::IMRMFeature* imrmfeature,
        #      const std::vector<std::string>& native_ids,
        #      const std::vector<double>& normalized_library_intensity,
        #      std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators,
        #      OpenSwath_Scores & scores) except + nogil 

        # void calculateLibraryScores(
        #      OpenSwath::IMRMFeature* imrmfeature,
        #      const std::vector<TransitionType> & transitions,
        #      const PeptideType& pep,
        #      const double normalized_feature_rt,
        #      OpenSwath_Scores & scores) except + nogil 

        # void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature, 
        #    const std::vector<TransitionType> & transitions,
        #    OpenSwath::SpectrumAccessPtr swath_map,
        #    OpenSwath::SpectrumAccessPtr ms1_map,
        #    OpenMS::DIAScoring & diascoring,
        #    const PeptideType& pep,
        #    OpenSwath_Scores & scores) except + nogil 

        void getNormalized_library_intensities_(libcpp_vector[LightTransition] transitions,
                                                libcpp_vector[double] normalized_library_intensity) except + nogil 

    cdef cppclass OpenSwath_Scores_Usage:

        OpenSwath_Scores_Usage() except + nogil 
        OpenSwath_Scores_Usage(OpenSwath_Scores_Usage) except + nogil  # wrap-ignore

        bool use_coelution_score_
        bool use_shape_score_
        bool use_rt_score_
        bool use_library_score_
        bool use_elution_model_score_
        bool use_intensity_score_
        bool use_total_xic_score_
        bool use_total_mi_score_
        bool use_nr_peaks_score_
        bool use_sn_score_
        bool use_mi_score_
        bool use_dia_scores_
        bool use_ms1_correlation
        bool use_ms1_fullscan
        bool use_ms1_mi
        bool use_uis_scores

    cdef cppclass OpenSwath_Scores:

        OpenSwath_Scores() except + nogil 
        OpenSwath_Scores(OpenSwath_Scores) except + nogil  # wrap-ignore

        double get_quick_lda_score(double library_corr_, double
                                   library_norm_manhattan_, double
                                   norm_rt_score_, double
                                   xcorr_coelution_score_, double
                                   xcorr_shape_score_, double log_sn_score_) except + nogil 

        double calculate_lda_prescore(OpenSwath_Scores scores) except + nogil 
        double calculate_swath_lda_prescore(OpenSwath_Scores scores) except + nogil 

        double elution_model_fit_score
        double library_corr
        double library_norm_manhattan
        double library_rootmeansquare
        double library_sangle
        double norm_rt_score

        double isotope_correlation
        double isotope_overlap
        double massdev_score
        double xcorr_coelution_score
        double xcorr_shape_score

        double yseries_score
        double bseries_score
        double log_sn_score

        double weighted_coelution_score
        double weighted_xcorr_shape
        double weighted_massdev_score
       
        double ms1_xcorr_coelution_score
        double ms1_xcorr_coelution_contrast_score
        double ms1_xcorr_coelution_combined_score
        double ms1_xcorr_shape_score
        double ms1_xcorr_shape_contrast_score
        double ms1_xcorr_shape_combined_score
        double ms1_ppm_score
        double ms1_isotope_correlation
        double ms1_isotope_overlap
        double ms1_mi_score
        double ms1_mi_contrast_score
        double ms1_mi_combined_score

        double library_manhattan
        double library_dotprod
        double intensity
        double total_xic
        double nr_peaks
        double sn_ratio
        double mi_score 
        double weighted_mi_score

        double rt_difference
        double normalized_experimental_rt
        double raw_rt_score

        double dotprod_score_dia
        double manhatt_score_dia
