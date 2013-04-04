// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
// **********************************************************************//
// CLASS Spec MERGER:
// merges 2 spectra which have been preprocessed
// by spec_merge-> common lc_peaks are marked!!!!
//

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_MS1FEATUREMERGER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_MS1FEATUREMERGER_H

namespace OpenMS
{

  class OPENMS_DLLAPI MS1FeatureMerger
  {


    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    //////
    // the common lc-ms spectrum,
    // created from teh overlap of A and B
    LCMS * lcmsMap;

    std::vector<int> idsToRemove;
    std::map<double, std::vector<SHFeature *> > mzClusters;


    ////////////////////////////////////////////////
    // declaration of the public members:

public:

    /*
    static double INTENSITY_APEX_THRESHOLD;
    static double MS1_PEAK_AREA_TR_RESOLUTION;

    static double INITIAL_TR_TOLERANCE;
    static double MS1_FEATURE_MERGING_TR_TOLERANCE;
    static double PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION;
    static double PPM_TOLERANCE_FOR_MZ_CLUSTERING;
    static bool MS1_FEATURE_CLUSTERING;
  */

    // class destructor
    ~MS1FeatureMerger();
    // class constructor
    MS1FeatureMerger(LCMS *);


    //////////////////////////////////////////////////
    // start the merging process
    void startFeatureMerging();
    // create a distribution of delta Tr for the splited features
    void createMZFeatureClusters();
    // process a vector of m/z features
    void processMZFeatureVector(std::vector<SHFeature *> *);
    // find to this feature the features which should be merged
    std::vector<SHFeature *>::iterator findFeaturesToMerge(SHFeature *, std::vector<SHFeature *>::iterator, std::vector<SHFeature *> *);
    // compare if a feature belongs to another feature
    bool compareMZFeatureBeloning(SHFeature *, SHFeature *);
    // merge the target to the search feature
    void mergeFeatures(SHFeature *, SHFeature *);
    // copmute new parameters for the merged MS1 feature
    void computeNewMS1FeatureParameters(SHFeature *);
    // computes the area of between 2 peaks:
    double computeDeltaArea(double, double, double, double);


    // this structure provides the function to compare
    // in the sorting algorithm:
    struct OPERATOR_FEATURE_TR
    {
      // provide the compare function for sort:
      bool operator()(const SHFeature A, const SHFeature B) const
      {
        // check if they have same mass
        return A.TR < B.TR;
      }

    };


    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

  };

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_MS1FEATUREMERGER_H
