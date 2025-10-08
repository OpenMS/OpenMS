// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/MATH/StatisticFunctions.h>

// #define DEBUG_MASSTRACES
// #include <assert.h>

bool SortDoubleDoublePairFirst(const std::pair<double, double>& left,
    const std::pair<double, double>& right);

namespace OpenMS
{
  using namespace std;
  using namespace OpenSwath;


  MasstraceCorrelator::MasstraceCorrelator()
  : DefaultParamHandler("MRMFeatureFinderScoring"),
    ProgressLogger()
  {   

    defaults_.setValue("sgolay_frame_length",15,"The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("sgolay_polynomial_order",3,"Order or the polynomial that is fitted.");
    defaults_.setValue("gauss_width",50,"Gaussian width.");

    // write defaults into Param object param_
    defaultsToParam_();
  }

  MasstraceCorrelator::~MasstraceCorrelator() = default;

  void MasstraceCorrelator::matchMassTraces_(
      const MasstracePointsType& hull_points1,
      const MasstracePointsType& hull_points2,
      std::vector<double>& vec1, std::vector<double>& vec2, double mindiff, double padEnds)
  {

    Size k=0,m=0;

    // If we do not pad the ends, we advance the longer array until the shorter one starts
    if (!padEnds)
    {
      while (k<hull_points1.size() && m<hull_points2.size() )
      {
        if (fabs(hull_points1[k].first - hull_points2[m].first) < mindiff)
        {
          break;
        }
        else if (hull_points1[k].first > hull_points2[m].first )
        {
            m++;
        }
        else if (hull_points1[k].first < hull_points2[m].first )
        {
            k++;
        }
      }
    
    }

    while (k<hull_points1.size() && m<hull_points2.size() )
    {
      if (fabs(hull_points1[k].first - hull_points2[m].first) < mindiff)
      {
        vec1.push_back(hull_points1[k].second);
        vec2.push_back(hull_points2[m].second);
        m++; k++;
      }
      else if (hull_points1[k].first > hull_points2[m].first )
      {
        //need to advance m, assume that vector 1 is zero
        vec1.push_back(0);
        vec2.push_back(hull_points2[m].second);
        m++;
      }
      else if (hull_points1[k].first < hull_points2[m].first )
      {
        //need to advance k, assume that vector 2 is zero
        vec1.push_back(hull_points1[k].second);
        vec2.push_back(0);
        k++;
      }
      else 
      {
        cout << "Error, cannot be here" << endl;
      }
    }

    // If we do not pad the ends, we can return now
    if (!padEnds) {return;}

    // If one vector is not at the end, we need to copy the rest and fill up with
    // zeros in the other.
    while (k<hull_points1.size())
    {
      vec1.push_back(hull_points1[k].second);
      vec2.push_back(0);
      k++;
    }

    while (m<hull_points2.size())
    {
      vec1.push_back(0);
      vec2.push_back(hull_points2[m].second);
      m++;
    }

  }

  void MasstraceCorrelator::scoreHullpoints(const MasstracePointsType& hull_points1, const MasstracePointsType& hull_points2,
        int& lag, double& lag_intensity, double& pearson_score, 
        const double min_corr, const int /* max_lag */, const double mindiff)
  {
    std::vector<double> vec1;
    std::vector<double> vec2;
    matchMassTraces_(hull_points1, hull_points2, vec1, vec2, mindiff);

    pearson_score = Math::pearsonCorrelationCoefficient(vec1.begin(), vec1.end(), vec2.begin(), vec2.end() );

    // If the correlation is below the minimum level, we can already return at this point
    if (pearson_score <= min_corr) 
    {
      return;
    }

    Scoring::XCorrArrayType xcorr_array = Scoring::normalizedCrossCorrelation(vec1, vec2, vec1.size(), 1);
    Scoring::XCorrArrayType::const_iterator pt = Scoring::xcorrArrayGetMaxPeak(xcorr_array);
    lag = pt->first;  // the lag / RT at the maximal Xcorr value =~ coelution score
    lag_intensity = pt->second; // the intensity at the maximal Xcorr value =~ shape score
  }

  void MasstraceCorrelator::createConsensusMapCache(const ConsensusMap& map, 
    std::vector< MasstracePointsType >& feature_points,
    std::vector< std::pair<double,double> >& max_intensities, 
    std::vector< double >& rt_cache)
  {

    startProgress(0, map.size(), "create consensus map cache");
    for (Size i = 0; i < map.size(); ++i)
    {
      setProgress(i);

      const ConsensusFeature::HandleSetType* f1_features = &map[i].getFeatures();

      // get the points into a vector of pairs (RT, intensity)
      MasstracePointsType f1_points; 
      for (ConsensusFeature::HandleSetType::iterator it = f1_features->begin(); it != f1_features->end(); ++it)
      {
        f1_points.push_back(std::make_pair(it->getRT(), it->getIntensity())); 
      }
      std::sort(f1_points.begin(), f1_points.end(), SortDoubleDoublePairFirst);
      feature_points.push_back(f1_points);

      // find maximum intensity and store it 
      double max_int = 0, max_mz =0;
      for (ConsensusFeature::HandleSetType::iterator it = f1_features->begin(); it != f1_features->end(); ++it)
      {
        if (it->getIntensity() > max_int)
        {
          max_int = it->getIntensity();
          max_mz  = it->getMZ();
        }
      }
      max_intensities.emplace_back(max_mz, max_int);
      rt_cache.push_back(map[i].getRT());
    }
    endProgress();
  
  }

  void MasstraceCorrelator::createPseudoSpectra(const ConsensusMap& map, 
      MSExperiment& pseudo_spectra,
      Size min_peak_nr, double min_correlation, 
      int max_lag, double max_rt_apex_difference)
  {

    // Parameters
    // double min_correlation = 0.7;
    // double max_lag = 1;
    // double max_rt_apex_difference = 3;
    // //double max_rt_apex_difference = 5000;

#ifdef DEBUG_MASSTRACES
    int opcounts = 0;
    int comparisons = 0;
    int nr_full_evals = 0;
    int nr_peaks_added = 0;
#endif

    Size j;
    double firstpoint, lastpoint, current_rt;
    double* rt_cache_ptr;

    int lag; double lag_intensity; double pearson_score;

    // cache datastructures
    std::vector< MasstracePointsType > feature_points; 
    std::vector< std::pair<double,double> > max_intensities; 
    std::vector< double > rt_cache;
    createConsensusMapCache(map, feature_points, max_intensities, rt_cache);

    std::map<int, int> used_already;
    // go through all consensus features in the map and use 
    startProgress(0, map.size(), "correlating masstraces ");
    for (Size i = 0; i < map.size(); ++i)
    {
      setProgress(i);

      if (used_already.find(i) != used_already.end()) 
      {
        continue;
      }
      used_already[i] = 0;

      // Prepare a new pseudo spectrum
      MSSpectrum spectrum;
      spectrum.getFloatDataArrays().clear();
      spectrum.getFloatDataArrays().resize(5);
      spectrum.getFloatDataArrays()[0].setName("RT_apex");
      spectrum.getFloatDataArrays()[1].setName("RT_diff");
      spectrum.getFloatDataArrays()[2].setName("lag");
      spectrum.getFloatDataArrays()[3].setName("pearson_score");
      spectrum.getFloatDataArrays()[4].setName("lag_intensity");
      spectrum.setRT(rt_cache[i]);
      spectrum.setMSLevel(2);

      // create the first peak of this spectrum == seed peak
      Peak1D peak;
      peak.setMZ(max_intensities[i].first);
      peak.setIntensity(max_intensities[i].second);
      spectrum.push_back(peak);

      // store the RT of the current feature and the first/last points of this feature
      firstpoint = feature_points[i].front().first; 
      lastpoint = feature_points[i].back().first;
      current_rt = rt_cache[i];
      rt_cache_ptr = &rt_cache[i];

      // go through all features with lower intensity in the map
      for (j = i+1; j < map.size(); ++j)
      {

        // If the center of this trace is outside the masstrace of the parent, ignore this pair.
        // If the difference between the rt_max of the two features is too large, ignore this pair.
        ++rt_cache_ptr;
        if ( fabs( (*rt_cache_ptr) - current_rt ) >  max_rt_apex_difference) {continue;}
        if ( (*rt_cache_ptr) < firstpoint || (*rt_cache_ptr) > lastpoint ) {continue;}

        // We score the two vectors against each other in terms of several properties / scores
        scoreHullpoints(feature_points[i], feature_points[j], lag, lag_intensity, pearson_score, min_correlation, max_lag);

#ifdef DEBUG_MASSTRACES
        cout << j << ". Checking mass trace at RT: "<<  map[j].getRT() << " m/z: " << map[j].getMZ()
          << " scores: [lag: " << lag << "] / [pearson: " << pearson_score << "]"<< endl;
        comparisons += 1;
        opcounts +=  feature_points[i].size() * feature_points[j].size();
        if (pearson_score > min_correlation)
        {
          nr_full_evals++;
        }
#endif

        // If all conditions are fulfilled, we add this feature as a peak. Note
        // that we need to check the pearson_score FIRST because the lag score is
        // only calculated if the pearson score is above the minimal value
        if (pearson_score > min_correlation && lag >= -max_lag && lag <= max_lag)
        {
          // mark this masstrace as used already, thus we cannot use it as a seed any more
          used_already[j] = 0;

#ifdef DEBUG_MASSTRACES
          nr_peaks_added++;
#endif

          Peak1D tmp_peak;
          tmp_peak.setMZ(max_intensities[j].first);
          tmp_peak.setIntensity(max_intensities[j].second);
          spectrum.push_back(tmp_peak);
          spectrum.getFloatDataArrays()[0].push_back(map[j].getRT());
          spectrum.getFloatDataArrays()[1].push_back(fabs(map[i].getRT() - map[j].getRT()));
          spectrum.getFloatDataArrays()[2].push_back(lag);
          spectrum.getFloatDataArrays()[3].push_back(pearson_score);
          spectrum.getFloatDataArrays()[4].push_back(lag_intensity);
        }
      }

      if (spectrum.size() > min_peak_nr)
      {
        pseudo_spectra.addSpectrum(spectrum);

#ifdef DEBUG_MASSTRACES
        cout << "Add spectrum " << i <<  " of size " << spectrum.size() << " at " << spectrum.getRT() << endl;
        cout << "===========================================================================" << endl;
#endif
      }

    }
    endProgress();

#ifdef DEBUG_MASSTRACES
    cout << "Nr operations " << opcounts << " / nr comparisons " << comparisons << " / full evaluations " << nr_full_evals << " :: nr spectra  " << pseudo_spectra.size() << endl;
    cout << "Nr peaks added " << nr_peaks_added << " out of total " << map.size() << endl;
#endif

  }


}

