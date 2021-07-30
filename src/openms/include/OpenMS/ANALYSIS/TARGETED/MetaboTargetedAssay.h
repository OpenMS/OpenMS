// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <map> //insert

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h> //SiriusMSFile
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h> //SiriusTargetDecoySpectra

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  /**
  @brief This class provides methods for the extraction of targeted assays for metabolomics.
  */

  class OPENMS_DLLAPI MetaboTargetedAssay
  {
    public:
    /**
    @brief MetaboTargetedAssay is able to store a precursor, metadata as well as compound information.
    */
    double precursor_int; ///< precursor intensity
    double transition_quality_score; ///< transitions quality score (not yet used)
    double precursor_mz; ///< precursor mass-to-charge
    double compound_rt; ///< compound retention time
    String molecular_formula; ///<< (putative) molecular formula
    int compound_file; ///< integer of file it belongs to in a list of files
    String compound_name; ///< compound name
    String compound_adduct; ///< compound adduct
    TargetedExperiment::Compound potential_cmp; ///< compound information stored in a TargetedExperiment
    std::vector<ReactionMonitoringTransition> potential_rmts; ///< vector of transitions belonging to the compound

    /**
    @brief CompoundTargetDecoyPair stores a pair of CompoundInfo and MSSpectrum (target, decoy)
    */
    class CompoundTargetDecoyPair
    {
    public:
      SiriusMSFile::CompoundInfo compound_info;
      SiriusFragmentAnnotation::SiriusTargetDecoySpectra target_decoy_spectra;

      CompoundTargetDecoyPair() = default;
      CompoundTargetDecoyPair(SiriusMSFile::CompoundInfo info, SiriusFragmentAnnotation::SiriusTargetDecoySpectra td_spectra) : compound_info(info), target_decoy_spectra(td_spectra) {}
    };

    /**
    @brief TargetDecoyGroup stores the mz, rt and file number in correspondence to the index of a MetaboTargetedAssay vector
    */

    class TargetDecoyGroup
    {
    public:
      int target_index = -1;
      int decoy_index = -1;
      double target_mz = 0.0;
      double target_rt = 0.0;
      double decoy_mz = 0.0;
      double decoy_rt = 0.0;
      int target_file_number = 0;
      int decoy_file_number = 0;
    };

    /**
    @brief Extract a vector of MetaboTargetedAssays without using fragment annotation

    @return Vector of MetaboTargetedAssay

    @param spectra Input of MSExperiment with spectra information
    @param feature_ms2_spectra_map FeatureMapping class with associated MS2 spectra
    @param precursor_rt_tol Retention time tolerance of the precursor
    @param precursor_mz_distance Max m/z distance of the precursor entries of two spectra to be merged
    @param cosine_sim_threshold Cosine similarity threshold for the usage of SpectraMerger
    @param transition_threshold Intensity threshold for MS2 peak used in MetaboTargetedAssay
    @param min_fragment_mz Minimum m/z a fragment ion has to have to be considered as a transition
    @param max_fragment_mz Maximum m/z a fragment ion has to have to be considered as a transition
    @param method_consensus_spectrum Boolean to use consensus spectrum method
    @param exclude_ms2_precursor Boolean to exclude MS2 precursor from MetaboTargetedAssay
    @param file_counter Count if multiple files are used.
    */
    static std::vector<MetaboTargetedAssay> extractMetaboTargetedAssay(const MSExperiment& spectra,
                                                                       const FeatureMapping::FeatureToMs2Indices& feature_ms2_index,
                                                                       const double& precursor_rt_tol,
                                                                       const double& precursor_mz_distance,
                                                                       const double& cosine_sim_threshold,
                                                                       const double& transition_threshold,
                                                                       const double& min_fragment_mz,
                                                                       const double& max_fragment_mz,
                                                                       const bool& method_consensus_spectrum,
                                                                       const bool& exclude_ms2_precursor,
                                                                       const unsigned int& file_counter);

    /**
    @brief Extract a vector of MetaboTargetedAssays using fragment annotation

    @return Vector of MetaboTargetedAssay

    @param v_cmp_spec Vector of CompoundInfo with associated fragment annotated MSspectrum
    @param transition_threshold Intensity threshold for MS2 peak used in MetaboTargetedAssay
    @param min_fragment_mz Minimum m/z a fragment ion has to have to be considered as a transition
    @param max_fragment_mz Maximum m/z a fragment ion has to have to be considered as a transition
    @param use_exact_mass Boolean if exact mass should be used as peak mass for annotated fragments
    @param exclude_ms2_precursor Boolean to exclude MS2 precursor from MetaboTargetedAssay
    @param file_counter Count if multiple files are used.
    */
    static std::vector<MetaboTargetedAssay> extractMetaboTargetedAssayFragmentAnnotation(const std::vector< CompoundTargetDecoyPair >& v_cmp_spec,
                                                                                         const double& transition_threshold,
                                                                                         const double& min_fragment_mz,
                                                                                         const double& max_fragment_mz,
                                                                                         const bool& use_exact_mass,
                                                                                         const bool& exclude_ms2_precursor,
                                                                                         const unsigned int& file_counter);

    /**
    @brief Pair compound information (SiriusMSFile) with the annotated target and decoy spectrum from SIRIUS/Passatutto based on the m_id (unique identifier composed of description_filepath_native_id_k introduced in the SiriusMSConverter)

    @return Vector of MetaboTargetedAssay::CompoundTargetDecoyPair

    @param v_cmpinfo Vector of SiriusMSFile::CompoundInfo
    @param annotated_spectra Vector of SiriusTargetDecoySpectra
    */
    static std::vector< MetaboTargetedAssay::CompoundTargetDecoyPair > pairCompoundWithAnnotatedSpectra(const std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                                                                                                        const std::vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra>& annotated_spectra);
    /**
    @brief Perform feature linking to build ambiguity groups based on the target and decoy position in the vector of MetaboTargetedAssays

    @return Map of pair (mz, rt) and vector of ambiguities for this mz,rt combination (MetaboTargetedAssay)

    @param v_mta Vector of MetaboTargetedAssay
    @param ar_mz_tol FeatureGroupingAlgorithmQT parameter distance_MZ:max_difference
    @param ar_rt_tol FeatureGroupingAlgorithmQT parameter distance_RT:max_difference
    @param ar_mz_tol_unit_res FeatureGroupingAlgorithmQT parameter distance_MZ_unit (ppm, Da)
    @param in_files_size Number of files which were processed in the vector of MetaboTargetedAssay (e.g. initially 5 different files in the vector<MetaboTargetedAssay>)
    */
    static std::unordered_map< UInt64, std::vector<MetaboTargetedAssay> > buildAmbiguityGroup(const std::vector<MetaboTargetedAssay>& v_mta,
                                                                                              const double& ar_mz_tol,
                                                                                              const double& ar_rt_tol,
                                                                                              const String& ar_mz_tol_unit_res, size_t in_files_size);

    /**
    @brief Resolve ambiguity groups based on occurrence in samples (e.g. at least in 20% of the samples) and if multiple possible identifications are reported within one ambiguity group use the one with the highest occurrence

    @return Map of pair (mz, rt) and vector of ambiguities for this mz,rt combination (MetaboTargetedAssay)

    @param total_occurrence_filter Value which has to be reached for the ambiguity group to be reported (e.g. in 20 % of the samples)
    @param in_files_size Number of files which were processed in the vector of MetaboTargetedAssay (e.g. initially 5 different files in the vector<MetaboTargetedAssay>)
    */
    static void resolveAmbiguityGroup(std::unordered_map< UInt64, std::vector<MetaboTargetedAssay> >& map_mta_filter,
                                      const double& total_occurrence_filter,
                                      size_t in_files_size);

  protected:

    /// Used to calculate the hard noise intensity threshold hard minimal threshold of min_int * noise_threshold_constant_
    static constexpr float noise_threshold_constant_ = 1.1;

    /**
    @brief Compare two peaks based on their intensity
    */
    static bool intensityLess_(Peak1D a, Peak1D b);

    /**
    @brief Gets charge from a singly charged adduct ([M+H]+/[M-H]-)
    */
    static int getChargeFromAdduct_(const String& adduct);

    /**
    @brief Filter one ambiguity group based on occurrence in samples (e.g. at least in 20% of the samples)

    @return Vector of MetaboTargetedAssay

    @param total_occurrence_filter Value which has to be reached for the ambiguity group to be reported (e.g. in 20 % of the samples)
    @param in_files_size Number of files which were processed in the vector of MetaboTargetedAssay (e.g. initially 5 different files in the vector<MetaboTargetedAssay>)
    */
    static void filterBasedOnTotalOccurrence_(std::vector<MetaboTargetedAssay>& mta, double total_occurrence_filter, size_t in_files_size);

    /**
    @brief Filter one ambiguity group with multiple possible identifications to use the one with the highest occurrence

    @return Vector of MetaboTargetedAssay
    */
    static void filterBasedOnMolFormAdductOccurrence_(std::vector<MetaboTargetedAssay>& mta);

    /**
    @brief Sort vector of MetaboTargetedAssay by precursor ion intensity
    */
    static void sortByPrecursorInt(std::vector<MetaboTargetedAssay>& vec_mta);
  };

} // namespace OpenMS
