// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief Feature trace in mass dimension for FLASHDeconv
  This class performs mass tracing on the deconvolved masses by FLASHDeconvAlgorithm
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

    /// Obtain and store information from deconvolved_spectrum (necessary information for mass tracing afterwards)
    void storeInformationFromDeconvolvedSpectrum(DeconvolvedSpectrum& deconvolved_spectrum);

    /**
       @brief Find mass features.
       @param averagine precalculated averagine for cosine calculation
       */
    std::vector<FLASHDeconvHelperStructs::MassFeature> findFeatures(const PrecalculatedAveragine& averagine);

  protected:
    void updateMembers_() override;

  private:
    /// cosine thresholds for scoring and filtering
    double min_isotope_cosine_;
    /// peak group information is stored in here for tracing
    std::map<double, std::map<double, PeakGroup>> peak_group_map_; // rt , mono mass, peakgroup
  };
} // namespace OpenMS