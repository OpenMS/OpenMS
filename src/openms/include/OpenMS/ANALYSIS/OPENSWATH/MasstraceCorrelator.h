// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MASSTRACECORRELATOR_H 
#define OPENMS_KERNEL_MASSTRACECORRELATOR_H 

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <fstream>

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left,
    const std::pair<double, double>& right);

namespace OpenMS
{


  /**
   @brief The MasstraceCorrelator tries to correlate individual masstraces found in mass spectrometric maps

   It does so using the normalized Cross-Correlation scoring of the OpenSWATH module.

   */
  class OPENMS_DLLAPI MasstraceCorrelator : 
    public DefaultParamHandler,
    public ProgressLogger
  {

  public:

    MasstraceCorrelator();

    ~MasstraceCorrelator();

    // a mass trace is a vector of pairs in (RT, Intensity)
    typedef std::vector<std::pair<double, double> > masstracePointsType;

    /** Compute pseudo-spectra from a set of (MS2) masstraces
     *
     * This function will take a set of masstraces (consensus map) as input and
     * produce a vector of pseudo spectra as output (pseudo_spectra result
     * vector).
     *
     * It basically makes an all-vs-all comparison of all masstraces against
     * each other and scores them on how similar they are in their elution
     * profiles.
     *
     * This assumes that the consensus feature is only from one (SWATH) map
     * This assumes that the consensus map is sorted by intensity
     *
     * (maybe) @todo use template peak type
     *
    */
    void createPseudoSpectra(ConsensusMap& map, MSExperiment<Peak1D>& pseudo_spectra,
        Size min_peak_nr, double min_correlation, int max_lag,
        double max_rt_apex_difference);

    /* Score two elution profiles against each other
     *
     * This function scores two elution profiles against each other:
     * First it creates 2 arrays that contain the corresponding intensities in RT-space
     * Then these arrays are scored using Cross-correlation scores and pearson coefficients.
     *
     * The pairs need to be sorted by the first entry (RT)
    */
    void scoreHullpoints(std::vector<std::pair<double, double> >& hull_points1,
        std::vector<std::pair<double, double> >& hull_points2, int& lag,
        double& lag_intensity, double& pearson_score, double min_corr,
        int max_lag, double mindiff = 0.1);

    /* Create a cache of the features in a consensus map
     *
     * This creates a cache of the input consensus map by creating the
     * following data structures:
     *  - a vector of mass traces (each mass trace is simply a vector of <RT,Intensity>
     *  - a vector of maximal intensities (max_rt, max_int)
     *  - a vector of retention times of the feature
     *
    */
    void createConsensusMapCache(const ConsensusMap& map,
        std::vector<masstracePointsType>& feature_points,
        std::vector<std::pair<double, double> >& max_intensities,
        std::vector<double>& rt_cache);

  private:

    /** @brief Match up two elution profile arrays
     *
     * This takes two 2D arrays of pairs (RT, int) and matches the entries to each
     * other on basis of RT, writing the corresponding intensities in the two
     * output vectors. If the RTs (pair.first are less than mindiff apart, the two
     * entries are considered to be equal, otherwise one is assumed to be zero).
     * This is useful for matching elution profiles that are not of the exact same
     * length and/or have missing values.
    */
    void match_elution_arrays(std::vector<std::pair<double, double> >& hull_points1,
        std::vector<std::pair<double, double> >& hull_points2,
        std::vector<double>& vec1, std::vector<double>& vec2, double mindiff,
        double padEnds = true);



  };
}

#endif
