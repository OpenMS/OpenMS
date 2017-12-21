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
// $Maintainer: Timo Sachsenberg $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSH_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSH_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnConfig.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

namespace OpenMS
{
  /**
      @brief The Superhirn FeatureFinderAlgorithm.

      The SuperHirn FeatureFinder algorithm is applied by calling "run" on this
      class, which in turn calls the FeatureFinderAlgorithmSHCtrl to execute
      the following algorithm:

       START by feeding the datavector into startScanParsing (FTPeakDetectController.cpp)
         For each scan
           1. Centroid it (new CentroidData instance), centroiding is done in
              CentroidData::calcCentroids
           2. Call add_scan_raw_data of ProcessData -> this also does the
              de-isotoping / feature finding in ProcessData::add_scan_raw_data
         3. Apply process_MS1_level_data_structure to the whole map
         4. Apply feature merging (MS1FeatureMerger)**, optionally
         5. Add to all LC MS/MS runs

         Step 2 in ProcessData::add_scan_raw_data works on centroided peaks of
          a single spectrum
           2.1 add to the background intensity controller
               BackgroundControl::addPeakMSScan which calculates intensity bins
           2.2 call "go" on Deisotoper (on single spectrum level) to "de-isotope" spectra *
           2.3 Converts them to objects of MSPeak type (single spectrum features)

         Step 3 works on an instance of ProcessData (clustering de-isotoped
          peaks from single spectra over RT) and applies the following steps:
           3.1 Extract elution peaks (call extract_elution_peaks of ProcessData)
           3.2 For all features, it creates a SuperHirn Feature (SHFeature)
           3.3 For all features, it computes the elution profile
               (FeatureLCProfile instance) and adds individual peaks to it

         Step 3.1 calls processIntensityMaps from BackgroundController

         * Deisotoper (Step 2.2)
           The Deisotoper works on single "peak groups" which is a set of
           peaks that has a maximal spacing of 1 + exps Da. These peak
           groups are produced by the CentroidData object
           [CentroidData::getNextPeakGroup] which internally holds a pointer
           to the current peak. It basically starts with the first peaks and
           adds peaks until the next peak is further away than 1+eps
           The Deisotoper then goes through the peak list, for each charge
           checks which peaks matches the current charge using
           IsotopicDist::getMatchingPeaks, creates an instance of DeconvPeak
           using this mono isotopic charge and then subtracts the current
           monoisotopic peak from the set of peaks using
           IsotopicDist::subtractMatchingPeaks (probably to account for
           overlapping isotopic patterns).
        ** Feature Merging Step 4 in MS1FeatureMerger::startFeatureMerging()
           which calls createMZFeatureClusters(). This method tries to check
           whether a feature is inside another feature using
           MS1FeatureMerger::compareMZFeatureBeloning which checks whether the ppm
           tolerance is below a certain level, the charge state is equal and
           whether both features have elution profiles.

      @ingroup FeatureFinder
  */
  class SUPERHIRN_DLLAPI FeatureFinderAlgorithmSH :
    public FeatureFinderAlgorithm,
    public FeatureFinderDefs
  {

public:
    typedef Peak1D PeakType;
    typedef FeatureFinderAlgorithm::MapType MapType; // MSExperiment
    typedef MapType::SpectrumType SpectrumType;

    FeatureFinderAlgorithmSH();

    unsigned int getNativeScanId(String native_id);

    void run() override;

    static FeatureFinderAlgorithm* create();

    static const String getProductName();

protected:
    MapType map_;
    using FeatureFinderAlgorithm::features_;

  };

}

#endif
