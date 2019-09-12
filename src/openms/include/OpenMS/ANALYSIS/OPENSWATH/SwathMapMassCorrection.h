// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
    ~SwathMapMassCorrection() = default;
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
     * @param swath_maps The raw swath maps from the current run
     *
     */
    void correctMZ(const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *>& transition_group_map,
                   std::vector< OpenSwath::SwathMap > & swath_maps, const OpenSwath::LightTargetedExperiment& targeted_exp);

    /**
     * @brief Correct the ion mobility values of a SWATH map based on the RT-normalization peptides
     *
     * This extracts the full spectra at the most likely elution time of the
     * calibrant peptides and fits a linear regression curve to correct for a
     * possible ion mobility (drift time) shift of the empirical drift time vs
     * the theoretically expected drift time.
     *
     * @param transition_group_map A MRMFeatureFinderScoring result map
     * @param swath_maps The raw swath maps from the current run
     * @param im_trafo The resulting map containing the transformation
     *
     */
    void correctIM(const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
                   const std::vector< OpenSwath::SwathMap > & swath_maps,
                   TransformationDescription& im_trafo,
                   const OpenSwath::LightTargetedExperiment& targeted_exp);

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

