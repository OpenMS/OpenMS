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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------
 
#pragma once

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
using namespace std;

namespace OpenMS
{
  /**
    description
  */
  class OPENMS_DLLAPI FIAMSDataProcessor
  {
public:
    /// Default constructor
    FIAMSDataProcessor();
    /// Default desctructor
    ~FIAMSDataProcessor();

    /// Copy constructor
    FIAMSDataProcessor(const FIAMSDataProcessor& cp);

    /// Assignment
    FIAMSDataProcessor& operator=(const FIAMSDataProcessor& fdp);

    // /// Process the given file
    // void process(std::string & filename, std::string & polarity, double resolution, double n_seconds);

    // /// Apply sliding bin smoothing
    // void applySlidingBin_(const MSSpectrum & input, MSSpectrum & output);

    /// Cut the spectra for time
    void cutForTime(const MSExperiment & experiment, vector<MSSpectrum> & output, float n_seconds);

    // /// Estimate noise for each peak
    // void trackNoise_(const MSSpectrum & input, MSSpectrum & output);

    // /// Perform accurate mass search
    // void runAccurateMassSearch_(const MSSpectrum & input, std::string & polarity, OpenMS::MzTab & output);
  };

} // namespace OpenMS

