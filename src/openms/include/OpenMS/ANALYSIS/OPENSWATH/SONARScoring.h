// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

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
    ~SONARScoring() override {}
    //@}
    //

    void computeSonarScores(OpenSwath::IMRMFeature* imrmfeature,
                            const std::vector<OpenSwath::LightTransition> & transitions,
                            const std::vector<OpenSwath::SwathMap>& swath_maps,
                            OpenSwath_Scores & scores) const;

private:

    void computeXCorr_(std::vector<std::vector<double> >& sonar_profiles,
                       double& xcorr_coelution_score, double& xcorr_shape_score) const;

    /// Copy constructor (algorithm class)
    SONARScoring(const SONARScoring& rhs);

    /// Assignment operator (algorithm class)
    SONARScoring& operator=(const SONARScoring& rhs);

    /// Synchronize members with param class
    void updateMembers_() override;

    double dia_extract_window_;
    bool dia_centroided_;
    bool dia_extraction_ppm_;

  };
}


