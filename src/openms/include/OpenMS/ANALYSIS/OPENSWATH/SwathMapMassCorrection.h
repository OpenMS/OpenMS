// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>

namespace OpenMS
{

  /**
   * @brief A class containing correction functions for Swath MS maps
   *
   * This class can use a set of pre-determined points in a Swath-MS map to
   * correct all maps according to the m/z shift found in those fixed points.
   *
   */
  class OPENMS_DLLAPI SwathMapMassCorrection :
    public DefaultParamHandler
  {

public:

    //@{
    /// Constructor
    SwathMapMassCorrection();

    /// Destructor
    ~SwathMapMassCorrection() override = default;
    //@}

    /// Synchronize members with param class
    void updateMembers_() override;

    /**
     * @brief Correct the m/z values of a SWATH map based on the RT-normalization peptides
     *
     * This extracts the full spectra at the most likely elution time of the
     * calibrant peptides and fits a regression curve to correct for a possible
     * mass shift of the empirical masses vs the theoretically expected masses.
     * Several types of regressions are available (see below corr_type parameter).
     *
     * The function will replace the pointers stored in swath_maps with a
     * transforming map that will contain corrected m/z values.
     *
     * @param transition_group_map A MRMFeatureFinderScoring result map
     * @param targeted_exp The corresponding spectral library (required for extraction coordinates)
     * @param swath_maps The raw swath maps from the current run, will be modified (replaced with a corrected version)
     * @param pasef Whether the data is PASEF data with possible overlapping m/z windows (with different ion mobility). In this case, the "best" SWATH window (with precursor cetntered around IM) is chosen.
     */
    void correctMZ(const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *>& transition_group_map,
                   const OpenSwath::LightTargetedExperiment & targeted_exp,
                   std::vector< OpenSwath::SwathMap > & swath_maps, const bool pasef);

    /**
     * @brief Correct the ion mobility values of a SWATH map based on the RT-normalization peptides
     *
     * This extracts the full spectra at the most likely elution time of the
     * calibrant peptides and fits a linear regression curve to correct for a
     * possible ion mobility (drift time) shift of the empirical drift time vs
     * the theoretically expected drift time. The resulting linear
     * transformation is stored using a TransformationDescription object.
     *
     * @param transition_group_map A MRMFeatureFinderScoring result map
     * @param swath_maps The raw swath maps from the current run
     * @param targeted_exp The corresponding spectral library (required for extraction coordinates)
     * @param pasef whether the data is PASEF data with possible overlapping m/z windows (with different ion mobility). In this case, the "best" SWATH window (with precursor cetntered around IM) is chosen.
     * @param im_trafo The resulting map containing the transformation
     */
    void correctIM(const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
                   const OpenSwath::LightTargetedExperiment & targeted_exp,
                   const std::vector< OpenSwath::SwathMap > & swath_maps,
                   const bool pasef,
                   TransformationDescription & im_trafo);

    /**
     * @brief Computes the SwathMaps for PASEF data in which windows can have the same m/z but differ by ion mobility
     *
     * For each precursor, the SwathMap is chosen based on library m/z and ion mobility.
     * If two or more SwathMaps isolate the same precursor the SwathMap in which the precursor is more centered across
     * ion mobility is chosen. Note that the single swath_map returned is in a vector so that this function is compatible with SONAR functions
     *
     * @param [IN] transition_group A MRMTransitionGroup for which the SwathMap is assigned to
     * @param [OUT] swath_maps A vector containing the a single entry, the swath map which the MRMFeature is assigned to
     */
    std::vector<OpenSwath::SwathMap> findSwathMapsPasef(const OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType& transition_group,
                                                         const std::vector< OpenSwath::SwathMap > & swath_maps);

  private:
    double mz_extraction_window_;
    bool mz_extraction_window_ppm_;
    bool ms1_im_;
    double im_extraction_window_;
    String mz_correction_function_;
    String im_correction_function_;
    String debug_im_file_;
    String debug_mz_file_;

  };
}

