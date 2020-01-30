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
#include <OpenMS/FORMAT/MzMLFile.h>

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
    FIAMSDataProcessor(
      String filename, 
      String dir_input, 
      String dir_output, 
      float resolution, 
      String polarity, 
      String db_mapping, 
      String db_struct, 
      String positive_adducts, 
      String negative_adducts, 
      bool store_progress=true,
      float min_mz=50, 
      float max_mz=1500, 
      float bin_step=20
    );

    /// Default desctructor
    ~FIAMSDataProcessor();

    /// Copy constructor
    FIAMSDataProcessor(const FIAMSDataProcessor& cp);

    /// Assignment
    FIAMSDataProcessor& operator=(const FIAMSDataProcessor& fdp);

    /// Process the given file
    void run(float n_seconds, OpenMS::MzTab & output);

    /// Cut the spectra for time
    void cutForTime(const MSExperiment & experiment, vector<MSSpectrum> & output, float n_seconds=6000);

    /// Merge spectra from different retention times into one
    MSSpectrum mergeAlongTime(const std::vector<OpenMS::MSSpectrum> & input);

    /// Pick peaks from merged spectrum and return as featureMap with the corresponding polarity
    MSSpectrum extractPeaks(const MSSpectrum & input);

    /// Convert a spectrum to a feature map with the corresponding polarity
    FeatureMap convertToFeatureMap(const MSSpectrum & input);

    /// Estimate noise for each peak
    MSSpectrum trackNoise(const MSSpectrum & input);

    /// Perform accurate mass search
    void runAccurateMassSearch(FeatureMap & input, OpenMS::MzTab & output);

    /// Get filename
    const String getFilename();

    /// Get input directory
    const String getInputDir();

    /// Get output directory
    const String getOutputDir();

    /// Get resolution
    const float getResolution();

    /// Get resolution
    const String getPolarity();

    /// Get minimum mass-to-charge
    const float getMinMZ();

    /// Get maximum mass-to-charge
    const float getMaxMZ();

    /// Get the sliding bin step
    const float getBinStep();

    /// Get the path to the db:mapping for passing to AccurateMassSearch
    const String getDBMapping();

    /// Get the path to the db:struct for passing to AccurateMassSearch
    const String getDBStruct();

    /// Get the path to the positive adducts for passing to AccurateMassSearch
    const String getPositiveAdducts();

    /// Get the path to the negative adducts for passing to AccurateMassSearch
    const String getNegativeAdducts();

    /// Get mass-to-charge ratios to base the sliding window upon
    const std::vector<float> getMZs();

    /// Get the sliding bin sizes
    const std::vector<float> getBinSizes();

private:
    void loadExperiment_();
    void storeSpectrum_(const MSSpectrum & input, String filename);

    String filename_;
    String dir_input_;
    String dir_output_;
    float resolution_;
    String polarity_;
    float min_mz_;
    float max_mz_;
    float bin_step_;
    String db_mapping_;
    String db_struct_;
    String positive_adducts_;
    String negative_adducts_;
    bool store_progress_;
    std::vector<float> mzs_; 
    std::vector<float> bin_sizes_;
    MSExperiment experiment_;
  };

} // namespace OpenMS

