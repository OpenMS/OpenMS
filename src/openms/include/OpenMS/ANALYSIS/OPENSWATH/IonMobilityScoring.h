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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// data access
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

// scoring
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <vector>

namespace OpenMS
{

  /** @brief A class that calls the ion mobility scoring routines
   *
   * Use this class to invoke the individual OpenSWATH ion mobility scoring routines.
   *
  */
  class OPENMS_DLLAPI IonMobilityScoring
  {
    typedef OpenSwath::LightCompound CompoundType;
    typedef OpenSwath::LightTransition TransitionType;

  public:

    /// Constructor
    IonMobilityScoring();

    /// Destructor
    ~IonMobilityScoring();

    /**
      @brief Performs scoring of the ion mobility dimension in MS2

      @param spectrum The DIA MS2 spectrum found at the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_lower Ion Mobility extraction start
      @param drift_upper Ion Mobility extraction end
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)

      @return Populates additional scores in the @p scores object

    */
    static void driftScoring(OpenSwath::SpectrumPtr spectrum,
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_lower,
                             const double drift_upper,
                             const double drift_target,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const bool use_spline,
                             const double drift_extra);

    /**
      @brief Performs scoring of the ion mobility dimension in MS1

      @param spectrum The DIA MS1 spectrum found at the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_lower Ion Mobility extraction start
      @param drift_upper Ion Mobility extraction end
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extra extraction to use for drift time (in percent)

      @return Populates additional scores in the @p scores object

    */
    static void driftScoringMS1(OpenSwath::SpectrumPtr spectrum,
                                const std::vector<TransitionType> & transitions,
                                OpenSwath_Scores & scores,
                                const double drift_lower,
                                const double drift_upper,
                                const double drift_target,
                                const double dia_extract_window_,
                                const bool dia_extraction_ppm_,
                                const bool use_spline,
                                const double drift_extra);

    /**
      @brief Performs scoring of the ion mobility dimension in MS1 and MS2 (contrast)

      @param spectrum The DIA MS2 spectrum found at the peak apex
      @param ms1spectrum The DIA MS1 spectrum found at the peak apex
      @param transitions The transitions used for scoring
      @param scores The output scores
      @param drift_lower Ion Mobility extraction start
      @param drift_upper Ion Mobility extraction end
      @param drift_target Ion Mobility extraction target
      @param dia_extraction_window_ m/z extraction width
      @param dia_extraction_ppm_ Whether m/z extraction width is in ppm
      @param use_spline Whether to use spline for fitting
      @param drift_extra Extra extraction to use for drift time (in percent)

      @return Populates additional scores in the @p scores object

    */
    static void driftScoringMS1Contrast(OpenSwath::SpectrumPtr spectrum, OpenSwath::SpectrumPtr ms1spectrum,
                                        const std::vector<TransitionType> & transitions,
                                        OpenSwath_Scores & scores,
                                        const double drift_lower,
                                        const double drift_upper,
                                        const double dia_extract_window_,
                                        const bool dia_extraction_ppm_,
                                        const double drift_extra);
  };
}

