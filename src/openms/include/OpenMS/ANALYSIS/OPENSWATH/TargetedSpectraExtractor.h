// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  /**
    @brief This class filters, annotates, picks, and scores spectra (e.g., taken
    from a DDA experiment) based on a target list.

    The input experiment is expected to be a standard .mzML file.
    The input target list is expected to be a standard TraML formatted .csv file.

    The class deals with filtering and annotating only those spectra that reflect
    the characteristics described by the target list.
    The filtering is done based on the transition_name, PrecursorMz and RetentionTime.
    The spectra are smoothed and peaks are picked for each spectrum.
    The spectra are then scored based on total TIC, FWHM and SNR.
    One spectrum is chosen for each of those transitions for which at least one
    valid spectrum was found and matched.

    The user can decide to use only extractSpectra(), otherwise run the methods
    in the following order:
    annotateSpectra()
    pickSpectrum() (called once for each annotated spectrum)
    scoreSpectra()
    selectSpectra()
  */
  class OPENMS_DLLAPI TargetedSpectraExtractor :
    public DefaultParamHandler
  {
public:
    TargetedSpectraExtractor();
    ~TargetedSpectraExtractor() override = default;

    /**
      Structure for a match against a spectral library
      TODO: Replace with MzTab once the final implementation is done
    */
    struct Match
    {
      Match() = default;
      Match(MSSpectrum a, double b) : spectrum(std::move(a)), score(b) {}
      MSSpectrum spectrum;
      double score = 0.0;
    };

    class OPENMS_DLLAPI Comparator
    {
    public:
      virtual ~Comparator() = default;
      virtual void generateScores(
        const MSSpectrum& spec,
        std::vector<std::pair<Size,double>>& scores,
        double min_score
      ) const = 0;

      virtual void init(
        const std::vector<MSSpectrum>& library,
        const std::map<String,DataValue>& options
      ) = 0;

      const std::vector<MSSpectrum>& getLibrary() const
      {
        return library_;
      }

    protected:
      std::vector<MSSpectrum> library_;
    };

    class OPENMS_DLLAPI BinnedSpectrumComparator : public Comparator
    {
    public:
      ~BinnedSpectrumComparator() override = default;
      void generateScores (
        const MSSpectrum& spec,
        std::vector<std::pair<Size,double>>& scores,
        double min_score
      ) const override
      {
        scores.clear();
        const BinnedSpectrum in_bs(spec, bin_size_, false, peak_spread_, bin_offset_);
        for (Size i = 0; i < bs_library_.size(); ++i)
        {
          const double cmp_score = cmp_bs_(in_bs, bs_library_[i]);
          if (cmp_score >= min_score)
          {
            scores.emplace_back(i, cmp_score);
          }
        }
      }

      void init(const std::vector<MSSpectrum>& library, const std::map<String,DataValue>& options) override;

    private:
      BinnedSpectralContrastAngle cmp_bs_;
      std::vector<BinnedSpectrum> bs_library_;
      double bin_size_ = 0.02; // Default for nominal mass: 1.0;
      UInt peak_spread_ = 0;
      double bin_offset_ = 0.0; // Default for nominal mass resolution: 0.4;
    };

    void getDefaultParameters(Param& params) const;

    /**
      @brief Filters and annotates those spectra that could potentially match the
      transitions of the target list.

      The spectra taken into account are those that fall within the precursor RT
      window and MZ tolerance set by the user through the parameters "rt_window"
      and "mz_tolerance". Default values are provided for both parameters.

      @warning The picked spectrum could be empty, meaning no peaks were found.

      @param[in] spectra The spectra to filter
      @param[in] targeted_exp The target list
      @param[out] annotated_spectra The spectra annotated with the related transition's name
      @param[out] features A FeatureMap which will contain informations about the name, precursor RT and precursor MZ of the matched transition
      @param[in] compute_features If false, `features` will be ignored
    */
    void annotateSpectra(
      const std::vector<MSSpectrum>& spectra,
      const TargetedExperiment& targeted_exp,
      std::vector<MSSpectrum>& annotated_spectra,
      FeatureMap& features,
      bool compute_features = true
    ) const;

    /**
      @brief Filters and annotates those spectra that could potentially match the
      transitions of the target list.

      The spectra taken into account are those that fall within the precursor RT
      window and MZ tolerance set by the user through the parameters "rt_window"
      and "mz_tolerance". Default values are provided for both parameters.

      @warning The picked spectrum could be empty, meaning no peaks were found.

      @param[in] spectra The spectra to filter
      @param[in] targeted_exp The target list
      @param[out] annotated_spectra The spectra annotated with the related transition's name
    */
    void annotateSpectra(
      const std::vector<MSSpectrum>& spectra,
      const TargetedExperiment& targeted_exp,
      std::vector<MSSpectrum>& annotated_spectra
    ) const;

    /**
      @brief Annotates the MS2 spectra with the likely MS1 feature that it was derived from

      Annotating based on MS1 feature results from AccurateMassSearch.
      In this case, the input will be e.g., const FeatureMap& ms1_features and the RTs and names (i.e., PeptideRef),
      defined in the FeatureMap.

      @param[in] spectra The spectra to filter
      @param[in] ms1_features The MS1 features
      @param[out] ms2_features The MS2 features
      @param[out] annotated_spectra The resulting annotated spectra
    */
    void annotateSpectra(
        const std::vector<MSSpectrum>& spectra,
        const FeatureMap& ms1_features,
        FeatureMap& ms2_features,
        std::vector<MSSpectrum>& annotated_spectra) const;
    
    /**
      @brief Search accurate masses and add identification (peptide hits) as features/sub-features

      @param[in] feat_map The feature map to search in
      @param[out] feat_map_output The output feature map, with peptide identifaction as sub features
      @param[in] add_unidentified_features Adds unidentified features to the feature map
    */
    void searchSpectrum(
        OpenMS::FeatureMap& feat_map,
        OpenMS::FeatureMap& feat_map_output,
        bool add_unidentified_features = false) const;

    /**
      @brief Picks a spectrum's peaks and saves them in picked_spectrum.

      The spectrum is first smoothed with a Gaussian filter (default) or using the
      Savitzky-Golay method. The parameter "use_gauss" handles this choice.
      The peak picking is executed with PeakPickerHiRes.
      Custom parameters provided with the prefix "GaussFilter:", "SavitzkyGolayFilter:"
      and "PeakPickerHiRes:" are taken into account.
      Peaks are then filtered by their heights and FWHM values.
      It is possible to set the peaks' limits through the parameters: "peak_height_min",
      "peak_height_max" and "fwhm_threshold".

      @throw Exception::IllegalArgument If `spectrum` is not sorted by position (mz).

      @param[in] spectrum The input spectrum
      @param[out] picked_spectrum A spectrum containing only the picked peaks
    */
    void pickSpectrum(const MSSpectrum& spectrum, MSSpectrum& picked_spectrum) const;

    /**
      @brief Assigns a score to the spectra given an input and saves them in scored_spectra.

      Also add the informations to the FeatureMap first constructed in annotateSpectra().
      The scores are based on total TIC, SNR and FWHM. It is possible to assign a
      weight to these parameters using: "tic_weight", "fwhm_weight" and "snr_weight".
      For each spectrum, the TIC and the SNR are computed on the entire spectrum.
      The FWHMs are computed only on picked peaks. Both SNR and FWHM are averaged values.
      The informations are added as FloatDataArray in scored_spectra and as MetaValue in features.

      @throw Exception::InvalidSize If `features` and `annotated_spectra` sizes don't match.

      @param[in] annotated_spectra The annotated spectra to score (for TIC and SNR)
      @param[in] picked_spectra The picked peaks found on each of the annotated spectra (for FWHM)
      @param[in,out] features The score informations are also added to this FeatureMap. Picked peaks' FWHMs are saved in features' subordinates.
      @param[out] scored_spectra The scored spectra. Basically a copy of annotated_spectra with the added score informations
      @param[in] compute_features If false, `features` will be ignored
    */
    void scoreSpectra(
      const std::vector<MSSpectrum>& annotated_spectra,
      const std::vector<MSSpectrum>& picked_spectra,
      FeatureMap& features,
      std::vector<MSSpectrum>& scored_spectra,
      bool compute_features = true
    ) const;

    /**
      @brief Assigns a score to the spectra given an input and saves them in scored_spectra.

      Also add the informations to the FeatureMap first constructed in annotateSpectra().
      The scores are based on total TIC, SNR and FWHM. It is possible to assign a
      weight to these parameters using: "tic_weight", "fwhm_weight" and "snr_weight".
      For each spectrum, the TIC and the SNR are computed on the entire spectrum.
      The FWHMs are computed only on picked peaks. Both SNR and FWHM are averaged values.

      @param[in] annotated_spectra The annotated spectra to score (for TIC and SNR)
      @param[in] picked_spectra The picked peaks found on each of the annotated spectra (for FWHM)
      @param[out] scored_spectra The scored spectra. Basically a copy of annotated_spectra with the added score informations
    */
    void scoreSpectra(
      const std::vector<MSSpectrum>& annotated_spectra,
      const std::vector<MSSpectrum>& picked_spectra,
      std::vector<MSSpectrum>& scored_spectra
    ) const;

    /**
      @brief The method selects the highest scoring spectrum for each possible
      annotation (i.e., transition name)

      @throw Exception::InvalidSize If `scored_spectra` and `features` sizes don't match.

      @param[in] scored_spectra Input annotated and scored spectra
      @param[in] features Input features
      @param[out] selected_spectra Output selected spectra
      @param[out] selected_features Output selected features
      @param[in] compute_features If false, `selected_features` will be ignored
    */
    void selectSpectra(
      const std::vector<MSSpectrum>& scored_spectra,
      const FeatureMap& features,
      std::vector<MSSpectrum>& selected_spectra,
      FeatureMap& selected_features,
      bool compute_features = true
    ) const;

    /**
      @brief The method selects the highest scoring spectrum for each possible
      annotation (i.e., transition name)

      @param[in] scored_spectra Input annotated and scored spectra
      @param[out] selected_spectra Output selected spectra
    */
    void selectSpectra(
      const std::vector<MSSpectrum>& scored_spectra,
      std::vector<MSSpectrum>& selected_spectra
    ) const;

    /**
      @brief Combines the functionalities given by all the other methods implemented
      in this class.

      The method expects an experiment and a target list in input,
      and constructs the extracted spectra and features.
      For each transition of the target list, the method tries to find its best
      spectrum match. A FeatureMap is also filled with informations about the
      extracted spectra.

      @param[in] experiment The input experiment
      @param[in] targeted_exp The target list
      @param[out] extracted_spectra The spectra related to the transitions
      @param[out] extracted_features The features related to the output spectra
      @param[in] compute_features If false, `extracted_features` will be ignored
    */
    void extractSpectra(
      const MSExperiment& experiment,
      const TargetedExperiment& targeted_exp,
      std::vector<MSSpectrum>& extracted_spectra,
      FeatureMap& extracted_features,
      bool compute_features = true
    ) const;

    /**
      @brief Combines the functionalities given by all the other methods implemented
      in this class.

      The method expects an experiment and a target list in input,
      and constructs the extracted spectra.
      For each transition of the target list, the method tries to find its best
      spectrum match.

      @param[in] experiment The input experiment
      @param[in] targeted_exp The target list
      @param[out] extracted_spectra The spectra related to the transitions
    */
    void extractSpectra(
      const MSExperiment& experiment,
      const TargetedExperiment& targeted_exp,
      std::vector<MSSpectrum>& extracted_spectra
    ) const;
    
    /**
      @brief Combines the functionalities given by all the other methods implemented
      in this class.

      The method expects an experiment and MS1 features in input,
      and constructs the extracted spectra and features.
      For each transition of the target list, the method tries to find its best
      spectrum match. A FeatureMap is also filled with informations about the
      extracted spectra.

      @param[in] experiment The input experiment
      @param[in] ms1_features The MS1 features map
      @param[out] extracted_spectra The spectra related to the transitions
    */
    void extractSpectra(
      const MSExperiment& experiment,
      const FeatureMap& ms1_features,
      std::vector<MSSpectrum>& extracted_spectra
    ) const;
    
    /**
      @brief Combines the functionalities given by all the other methods implemented
      in this class.

      The method expects an experiment and MS1 features in input,
      and constructs the extracted spectra and features.
      For each transition of the target list, the method tries to find its best
      spectrum match. A FeatureMap is also filled with informations about the
      extracted spectra.

      @param[in] experiment The input experiment
      @param[in] ms1_features The MS1 features map
      @param[out] extracted_spectra The spectra related to the transitions
      @param[out] extracted_features The features related to the output spectra
    */
    void extractSpectra(
      const MSExperiment& experiment,
      const FeatureMap& ms1_features,
      std::vector<MSSpectrum>& extracted_spectra,
      FeatureMap& extracted_features
    ) const; 
    
    /**
      @brief Searches the spectral library for the top scoring candidates that
      match the input spectrum.

      @param[in] input_spectrum The input spectrum for which a match is desired
      @param[in] cmp The comparator object containing the library and the logic for matching
      @param[out] matches A vector of `Match`es, containing the matched spectra and their scores
    */
    void matchSpectrum(
      const MSSpectrum& input_spectrum,
      const Comparator& cmp,
      std::vector<Match>& matches
    ) const;

    /**
      @brief Compares a list of spectra against a spectral library and updates
      the related features.

      The metavalues added to each `Feature` within the `FeatureMap` are:
      - spectral_library_name The name of the match's spectrum found in the library
      - spectral_library_score The match score [0-1]
      - spectral_library_comments The comments for the match's spectrum

      If a match for a given input spectrum is not found, the metavalues will be
      assigned a default value:
      - spectral_library_name and spectral_library_comments: an empty string
      - spectral_library_score: a value of 0.0

      @note The input `spectra` (and related `features`) are assumed to be the
      result of `extractSpectra()`, meaning they went (at least) through the process
      of peak picking.

      @param[in] spectra The input spectra
      @param[in] cmp The `Comparator` object containing the spectral library
      @param[in,out] features The `FeatureMap` to be updated with matching info
    */
    void targetedMatching(
      const std::vector<MSSpectrum>& spectra,
      const Comparator& cmp,
      FeatureMap& features
    );

    /**
      @brief Compares a list of spectra against a spectral library and creates
      a `FeatureMap` with the relevant information.

      The metavalues added to each `Feature` within the `FeatureMap` are:
      - spectral_library_name The name of the match's spectrum found in the library
      - spectral_library_score The match score [0-1]
      - spectral_library_comments The comments for the match's spectrum

      If a match for a given input spectrum is not found, the metavalues will be
      assigned a default value:
      - spectral_library_name and spectral_library_comments: an empty string
      - spectral_library_score: a value of 0.0

      @note The input `spectra` (and related `features`) are assumed to be unprocessed,
      therefore undergoing a process of peak picking during execution of this method.

      @param[in] spectra The input spectra
      @param[in] cmp The `Comparator` object containing the spectral library
      @param[out] features The `FeatureMap` to be filled with matching info
    */
    void untargetedMatching(
      const std::vector<MSSpectrum>& spectra,
      const Comparator& cmp,
      FeatureMap& features
    );

    /**
      @brief compute transitions list from MS1 and the associated MS2 features

      @param[in] ms1_features the MS1 features
      @param[in] ms2_features the MS2 features
      @param[out] t_exp the targeted experiment, containing the transitions
    */
    void constructTransitionsList(const OpenMS::FeatureMap& ms1_features, const OpenMS::FeatureMap& ms2_features, TargetedExperiment& t_exp) const;

    /**
      @brief store spectra in MSP format

      @param[in] filename the filename of the file to write
      @param[in] experiment the experiment to store
    */
    void storeSpectraMSP(const String& filename, MSExperiment& experiment) const;
    
    /**
      @brief organize into a map by combining features and subordinates with the same `identifier`

      @param[in] fmap_input input features map
      @param[in] fmap_output output features map
    */
    void mergeFeatures(const OpenMS::FeatureMap& fmap_input, OpenMS::FeatureMap& fmap_output) const;

protected:
    /// Overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
    void updateMembers_() override;

    /// Deisotope MS2 spectra
    void deisotopeMS2Spectra_(MSExperiment& experiment) const;

    /// Remove peaks form MS2 which are at a higher mz than the precursor + 10 ppm
    void removeMS2SpectraPeaks_(MSExperiment& experiment) const;

    /// organize into a map by combining features and subordinates with the same `identifier`
    void organizeMapWithSameIdentifier(const OpenMS::FeatureMap& fmap_input, std::map<OpenMS::String, std::vector<OpenMS::Feature>>& fmapmap) const;

  private:
    /**
      @brief Combines the functionalities given by all the other methods implemented
      in this class.

      The method expects an experiment and MS1 features in input,
      and constructs the extracted spectra and features.
      For each transition of the target list, the method tries to find its best
      spectrum match. A FeatureMap is also filled with informations about the
      extracted spectra.

      @param[in] experiment The input experiment
      @param[in] ms1_features The MS1 features map
      @param[out] extracted_spectra The spectra related to the transitions
      @param[out] extracted_features The features related to the output spectra
      @param[in] compute_features If false, `extracted_features` will be ignored
    */
    void extractSpectra(
      const MSExperiment& experiment,
      const FeatureMap& ms1_features,
      std::vector<MSSpectrum>& extracted_spectra,
      FeatureMap& extracted_features,
      const bool compute_features
    ) const;

  private:
    /**
      Unit to use for mz_tolerance_ and fwhm_threshold_: true for Da, false for ppm.
    */
    bool mz_unit_is_Da_;

    /**
      Precursor Retention Time window used during the annotation phase.
      For each transition in the target list, annotateSpectra() looks for
      the first spectrum whose RT time falls within the RT Window, whose
      left and right limits are computed at each analyzed spectrum.
      Also the spectrum's precursor MZ is checked against the transition MZ.
    */
    double rt_window_;

    /**
      Precursor MZ tolerance used during the annotation phase.
      For each transition in the target list, annotateSpectra() looks for
      the first spectrum whose precursor MZ is close enough (+-mz_tolerance_)
      to the transition's MZ.
      Also the spectrum's precursor RT is checked against the transition RT.
    */
    double mz_tolerance_;

    /**
      Used in pickSpectrum(), a peak's intensity needs to be >= peak_height_min_
      for it to be picked.
    */
    double peak_height_min_;

    /**
      Used in pickSpectrum(), a peak's intensity needs to be <= peak_height_max_
      for it to be picked.
    */
    double peak_height_max_;

    /**
      Used in pickSpectrum(), a peak's FWHM needs to be >= fwhm_threshold_
      for it to be picked.
    */
    double fwhm_threshold_;

    double tic_weight_; /**< Total TIC's weight when computing a spectrum's score */
    double fwhm_weight_; /**< FWHM's weight when computing a spectrum's score */
    double snr_weight_; /**< SNR's weight when computing a spectrum's score */

    /**
      Used in selectSpectra(), after the spectra have been assigned a score.
      Remained transitions will have at least one spectrum assigned.
      Each spectrum needs to have a score >= min_select_score_ to be valid,
      otherwise it gets filtered out.
    */
    double min_select_score_;

    /**
      Used in pickSpectrum(), it selects which filtering method is used during
      the smoothing phase.
      By default the Gauss filter is selected. Set to false for the Savitzky-Golay method.
    */
    bool use_gauss_;

    /**
      The number of matches to output from `matchSpectrum()`.
      These will be the matches of highest scores, sorted in descending order.
    */
    Size top_matches_to_report_;

    /// Minimum score for a match to be considered valid in `matchSpectrum()`.
    double min_match_score_;

    double min_fragment_mz_;
    double max_fragment_mz_;

    double relative_allowable_product_mass_;

    bool deisotoping_use_deisotoper_;
    double deisotoping_fragment_tolerance_;
    std::string deisotoping_fragment_unit_;
    int deisotoping_min_charge_;
    int deisotoping_max_charge_;
    int deisotoping_min_isopeaks_;
    int deisotoping_max_isopeaks_;
    bool deisotoping_keep_only_deisotoped_;
    bool deisotoping_annotate_charge_;

    double max_precursor_mass_threashold_;
  };
}
