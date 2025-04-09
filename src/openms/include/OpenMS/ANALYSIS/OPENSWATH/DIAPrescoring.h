// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

namespace OpenMS
{
  /**
    @brief Scoring of an spectrum given library intensities of a transition group.

    In DIA (data independent acquisition) / SWATH analysis, at each
    chromatographic point a full MS2 spectrum is recorded. This class allows to
    compute a number of scores based on the full MS2 spectrum available. The scores are the following:

    See also class DIAScoring.

    Simulate theoretical spectrum from library intensities of transition group
    and compute manhattan distance and dotprod score between spectrum intensities
    and simulated spectrum.
  */

  class OPENMS_DLLAPI DiaPrescore :
    public DefaultParamHandler
  {
    double dia_extract_window_; //done
    int nr_isotopes_;
    int nr_charges_;
public:

    DiaPrescore();

    DiaPrescore(double dia_extract_window, int nr_isotopes = 4, int nr_charges = 4);

    void defineDefaults();

    void updateMembers_() override;

    /**
      @brief Score a spectrum given a transition group.

      Simulate theoretical spectrum from library intensities of transition group
      and compute manhattan distance and dotprod score between spectrum intensities
      and simulated spectrum.
    */
    void score(const SpectrumSequence& spec,
               const std::vector<OpenSwath::LightTransition>& lt,
               const RangeMobility& im_range,
               double& dotprod,
               double& manhattan) const;

    /**
      @brief Compute manhattan and dotprod score for all spectra which can be accessed by
      the SpectrumAccessPtr for all transitions groups in the LightTargetedExperiment.
    */
    void operator()(const OpenSwath::SpectrumAccessPtr& swath_ptr,
                    OpenSwath::LightTargetedExperiment& transition_exp_used, const RangeMobility& range_im,
                    OpenSwath::IDataFrameWriter* ivw) const;
  };


}

