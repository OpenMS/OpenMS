// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_TARGETEDSPECTRAEXTRACTOR_H
#define OPENMS_ANALYSIS_OPENSWATH_TARGETEDSPECTRAEXTRACTOR_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

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
    virtual ~TargetedSpectraExtractor();

    void getDefaultParameters(Param& params);

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
    */
    void annotateSpectra(
      const std::vector<MSSpectrum>& spectra,
      const TargetedExperiment& targeted_exp,
      std::vector<MSSpectrum>& annotated_spectra,
      FeatureMap& features
    );

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

      @param[in] spectrum The input spectrum
      @param[out] picked_spectrum A spectrum containing only the picked peaks
    */
    void pickSpectrum(const MSSpectrum& spectrum, MSSpectrum& picked_spectrum);

    /**
      @brief Assigns a score to the spectra given an input and saves them in scored_spectra.

      Also add the informations to the FeatureMap first constructed in annotateSpectra().
      The scores are based on total TIC, SNR and FWHM. It is possible to assign a
      weight to these parameters using: "tic_weight", "fwhm_weight" and "snr_weight".
      For each spectrum, the TIC and the SNR are computed on the entire spectrum.
      The FWHMs are computed only on picked peaks. Both SNR and FWHM are averaged values.
      The informations are added as FloatDataArray in scored_spectra and as MetaValue in features.

      @param[in] annotated_spectra The annotated spectra to score (for TIC and SNR)
      @param[in] picked_spectra The picked peaks found on each of the annotated spectra (for FWHM)
      @param[in,out] features The score informations are also added to this FeatureMap. Picked peaks' FWHMs are saved in features' subordinates.
      @param[out] scored_spectra The scored spectra. Basically a copy of annotated_spectra with the added score informations
    */
    void scoreSpectra(
      const std::vector<MSSpectrum>& annotated_spectra,
      const std::vector<MSSpectrum>& picked_spectra,
      FeatureMap& features,
      std::vector<MSSpectrum>& scored_spectra
    );

    /**
      @brief The method selects the highest scoring spectrum for each possible
      annotation (i.e., transition name)

      @param[in] scored_spectra Input annotated and scored spectra
      @param[in] features Input features
      @param[out] selected_spectra Output selected spectra
      @param[out] selected_features Output selected features
    */
    void selectSpectra(
      const std::vector<MSSpectrum>& scored_spectra,
      const FeatureMap& features,
      std::vector<MSSpectrum>& selected_spectra,
      FeatureMap& selected_features
    );

    /**
      @brief The method selects the highest scoring spectrum for each possible
      annotation (i.e., transition name)

      @param[in] scored_spectra Input annotated and scored spectra
      @param[out] selected_spectra Output selected spectra
    */
    void selectSpectra(
      const std::vector<MSSpectrum>& scored_spectra,
      std::vector<MSSpectrum>& selected_spectra
    );

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
    */
    void extractSpectra(
      const MSExperiment& experiment,
      const TargetedExperiment& targeted_exp,
      std::vector<MSSpectrum>& extracted_spectra,
      FeatureMap& extracted_features
    );

protected:
    /// Overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed
    void updateMembers_();

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
      Also the spectrum's percursor MZ is checked against the transition MZ.
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
      Each spectrum needs to have a score >= min_score_ to be valid,
      otherwise it gets filtered out.
    */
    double min_score_;

    /**
      Used in pickSpectrum(), it selects which filtering method is used during
      the smoothing phase.
      By default the Gauss filter is selected. Set to false for the Savitzky-Golay method.
    */
    bool use_gauss_;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_TARGETEDSPECTRAEXTRACTOR_H
