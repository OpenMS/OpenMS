//--------------------------------------------------------------------------
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