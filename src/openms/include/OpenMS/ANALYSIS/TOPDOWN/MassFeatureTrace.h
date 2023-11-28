// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief Feature trace in mass dimension for FLASHDeconv
  This class performs mass tracing on the deconvolved masses by SpectralDeconvolution
  In other words, per spectrum deconvolved masses are converted into deconvolved features
  Currently only works for MS1 spectra. (Top-down DIA is not yet used much).
  Every time an MS1 spectrum is deconvolved, the relevant information is stored in this class.
  Tracing is performed at the end of FLASHDeconv run.
  This class also comes with tsv, TopFD, ProMex format output functions.
  @ingroup Topdown
  */

  class OPENMS_DLLAPI MassFeatureTrace : public DefaultParamHandler
  {
  public:
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// constructor
    MassFeatureTrace();

    /// destructor
    ~MassFeatureTrace() override = default;

    /// copy constructor
    MassFeatureTrace(const MassFeatureTrace&) = default;

    /// move constructor
    MassFeatureTrace(MassFeatureTrace&& other) = default;

    /// assignment operator
    MassFeatureTrace& operator=(const MassFeatureTrace& fd) = default;
    MassFeatureTrace& operator=(MassFeatureTrace&& fd) = default;

    /**
       @brief Find mass features.
       @param averagine precalculated averagine for cosine calculation
       @param ms_level ms level to process
       @param is_decoy if set, only process decoy spectra. otherwise only target spectra
       */
    std::vector<FLASHDeconvHelperStructs::MassFeature> findFeaturesAndUpdateQscore2D(const PrecalculatedAveragine& averagine, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level = 1,
                                                                                     bool is_decoy = false);

  protected:
    void updateMembers_() override;

  private:
    /// peak group information is stored in here for tracing
    std::map<double, std::map<double, PeakGroup>> peak_group_map_; // rt , mono mass, peakgroup
  };
} // namespace OpenMS