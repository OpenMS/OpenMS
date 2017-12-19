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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MASSTRACECORRELATOR_H 
#define OPENMS_ANALYSIS_OPENSWATH_MASSTRACECORRELATOR_H 

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
   @brief Correlates individual masstraces found in mass spectrometric maps

   The MasstraceCorrelator offers several functions to correlate individual
   mass traces using the normalized Cross-Correlation and pearson scoring of
   the OpenSWATH module.

   */
  class OPENMS_DLLAPI MasstraceCorrelator : 
    public DefaultParamHandler,
    public ProgressLogger
  {

  public:

    MasstraceCorrelator();

    ~MasstraceCorrelator() override;

    // a mass trace is a vector of pairs in (RT, Intensity)
    typedef std::vector<std::pair<double, double> > MasstracePointsType;

    /** Compute pseudo-spectra from a set of (MS2) masstraces
     *
     * This function will take a set of masstraces (consensus map) as input and
     * produce a vector of pseudo spectra as output (pseudo_spectra result
     * vector).
     *
     * It basically makes an all-vs-all comparison of all masstraces against
     * each other and scores them on how similar they are in their mass traces.
     *
     * This assumes that the consensus feature is only from one (SWATH) map
     * This assumes that the consensus map is sorted by intensity
     *
    */
    void createPseudoSpectra(ConsensusMap& map, MSExperiment& pseudo_spectra,
        Size min_peak_nr, double min_correlation, int max_lag,
        double max_rt_apex_difference);

    /* Score two mass traces against each other
     *
     * This function scores two mass traces (vector of <RT,Intensity>) against each other:
     *
     *  - The algorithm first creates 2 arrays that contain matched intensities
     *    in RT-space (accounting for missing data points and unequal length)
     *  - Next, these arrays are scored using cross-correlation scores and
     *    pearson coefficients.
     *
     * @note The pairs need to be sorted by the first entry (RT)
     *
     * @param hull_points1 The first input masstrace
     * @param hull_points2 The second input masstrace
     * @param lag The computed lag (output coelution score) 
     * @param lag_intensity The computed intensity at the lag (output shape score) 
     * @param pearson_score The computed pearson score (output)
     * @param min_corr Minimal correlation needed to proceed computing the cross-correlations
     * @param max_lag Currently unused
     * @param mindiff Minimal differences for matching up the two mass traces
     *
    */
    void scoreHullpoints(const MasstracePointsType& hull_points1, const MasstracePointsType& hull_points2,
        int& lag, double& lag_intensity, double& pearson_score, 
        const double min_corr, const int max_lag, const double mindiff = 0.1);

    /* Create a cache of the features in a consensus map
     *
     * This creates a cache of the input consensus map by creating the
     * following data structures:
     *  - a vector of mass traces (each mass trace is simply a vector of <RT,Intensity>
     *  - a vector of maximal intensities (max_rt, max_int)
     *  - a vector of retention times of the feature
     *
     * @param map The input consensus map
     * @param feature_points The list of all mass traces
     * @param max_intensities The list of maximal intensities 
     * @param rt_cache The list of retention times of all features
    */
    void createConsensusMapCache(const ConsensusMap& map,
        std::vector<MasstracePointsType>& feature_points,
        std::vector<std::pair<double, double> >& max_intensities,
        std::vector<double>& rt_cache);

  protected:

    /** @brief Match up two mass traces with potentially missing values
     *
     * To compute correlations on masstraces, they need to have the same length
     * and matching points. This function matches two masstraces by RT and
     * identifies points that are the same in retention time (see mindiff
     * parameter) and matches them. If no match is found, a missing value is
     * assumed and they are filled with zeros. Thus, if the two retention times
     * are less than mindiff apart, the two entries are considered to be equal,
     * otherwise one is assumed to be zero).
     *
     * This is useful for matching mass traces that are not of the exact same
     * length and/or have missing values.
     *
     * @param hull_points1 The first input mass trace
     * @param hull_points2 The second input mass trace
     * @param vec1 The intensities of the first mass trace with matched-up points
     * @param vec2 The intensities of the second mass trace with matched-up points
     * @param mindiff The minimal difference in RT for points to match up 
     * @param padEnds Whether to pad ends with zeros
     *
    */
    void matchMassTraces_(const MasstracePointsType& hull_points1, const MasstracePointsType& hull_points2,
        std::vector<double>& vec1, std::vector<double>& vec2, double mindiff,
        double padEnds = true);
  };
}

#endif
