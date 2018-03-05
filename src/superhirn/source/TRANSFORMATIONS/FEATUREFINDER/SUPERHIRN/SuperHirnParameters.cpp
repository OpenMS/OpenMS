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
// $Maintainer: Peter Kunszt $
// $Authors: Lukas Mueller, Markus Mueller Peter Kunszt $
// --------------------------------------------------------------------------
//

#include <map>
#include <cstdio>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

  SuperHirnParameters * SuperHirnParameters::instance_ = nullptr;
  bool SuperHirnParameters::haveInstance_ = false;

  SuperHirnParameters::SuperHirnParameters()
  {
    // Background Intensity static vars
    backgroundIntensityBinsTR_ = 2.0;
    backgroundIntensityBinsMZ_ = 50.0;
    backgroundIntensityBinsIntens_ = 50;
    backgroundIntensityBinsMinBinCount_ = 1;

    // these will be overwritten at config time
    minTR_ = 0;
    maxTR_ = 0;
    minFeatureMZ_ = 0;
    maxFeatureMZ_ = 0;
    minFeatureChrg_ = 0;
    maxFeatureChrg_ = 0;


    // minimal intensity level:  NEVER USED
    intensityThreshold_ = 0;

    // m/z tolerance value: NEVER CONFIGURED
    toleranceMZ_ = 10;

    // max_distance from next elution peak member in min.:
    maxInterScanRetentionTimeDistance_ = 0;

    // define minimal number of members in LC elution peaks cluster
    minNbClusterMembers_ = 0;

    /*
     // to track detected monoistopic mass for debugging:
     monoIsoDebugging_ = false;
     debugMonoIsoMassMin_ = 318.00;
     debugMonoIsoMassMax_ = 319.00;
     */

    // if data are in centroid form or not:
    centroidDataModus_ = false;

    massTolPpm_ = 10.0;
    massTolDa_ = 0.01;
    minIntensity_ = 0.0;
    intensityFloor_ = 1.0;

    peptideProbabilityThreshold_ = 0.9;     // this is hardcoded. it is never configured.
    storeAllLowProbabilityMS2Scans_ = false;     // this is hardcoded.


    createFeatureElutionProfiles_ = false;
    /*
    elutionPeakDebugging_ = false;
    elutionPeakMassMin_ = -1;
    elutionPeakMassMax_ = -2;

    // if this option is on, then construct fake features for available MS2 features
    featureFakeInsertionBasedOnMS2Feature_ = true;
    */

    lowIntensityMSSignalThreshold_ = 1.0;     // never configured, but used
    initIsotopeDist_ = false;
  }

}
