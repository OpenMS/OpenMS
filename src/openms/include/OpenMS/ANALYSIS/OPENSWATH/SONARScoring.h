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

#ifndef OPENMS_ANALYSIS_OPENSWATH_SONARSCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_SONARSCORING_H

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>

namespace OpenMS
{
  /**
    @brief Scoring of an spectrum using SONAR data

    @htmlinclude OpenMS_SONARScoring.parameters

  */
  class OPENMS_DLLAPI SONARScoring :
    public DefaultParamHandler
  {
public:

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    SONARScoring();

    /// Destructor
    virtual ~SONARScoring() {}
    //@}
    //

    void computeSonarScores(OpenSwath::IMRMFeature* imrmfeature,
                      const std::vector<OpenSwath::LightTransition> & transitions,
                      std::vector<OpenSwath::SwathMap>& swath_maps,
                      OpenSwath_Scores & scores);

private:

    void computeXCorr_(std::vector<std::vector<double> >& sonar_profiles,
                       double& xcorr_coelution_score, double& xcorr_shape_score);

    /// Copy constructor (algorithm class)
    SONARScoring(const SONARScoring& rhs);

    /// Assignment operator (algorithm class)
    SONARScoring& operator=(const SONARScoring& rhs);

    /// Synchronize members with param class
    void updateMembers_();

    double dia_extract_window_;
    bool dia_centroided_;
    bool dia_extraction_ppm_;

  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_SONARSCORING_H

