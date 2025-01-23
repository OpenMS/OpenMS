// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>

#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>

#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  LabeledPairFinder::LabeledPairFinder() :
    BaseGroupFinder()
  {
    setName("LabeledPairFinder");
    defaults_.setValue("rt_estimate", "true", "If 'true' the optimal RT pair distance and deviation are estimated by "
                                              "fitting a gaussian distribution to the histogram of pair distance. "
                                              "Note that this works only datasets with a significant amount of pairs! "
                                              "If 'false' the parameters 'rt_pair_dist', 'rt_dev_low' "
                                              "and 'rt_dev_high' define the optimal distance.");
    defaults_.setValidStrings("rt_estimate", {"true","false"});
    defaults_.setValue("rt_pair_dist", -20.0, "optimal pair distance in RT [sec] from light to heavy feature");
    defaults_.setValue("rt_dev_low", 15.0, "maximum allowed deviation below optimal retention time distance");
    defaults_.setMinFloat("rt_dev_low", 0.0);
    defaults_.setValue("rt_dev_high", 15.0, "maximum allowed deviation above optimal retention time distance");
    defaults_.setMinFloat("rt_dev_high", 0.0);

    defaults_.setValue("mz_pair_dists", ListUtils::create<double>("4.0"), "optimal pair distances in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)");
    defaults_.setValue("mz_dev", 0.05, "maximum allowed deviation from optimal m/z distance\n");
    defaults_.setMinFloat("mz_dev", 0.0);
    defaults_.setValue("mrm", "false", "this option should be used if the features correspond mrm chromatograms (additionally the precursor is taken into account)", {"advanced"});
    defaults_.setValidStrings("mrm", {"true","false"});

    defaultsToParam_();
  }

  void LabeledPairFinder::run(const vector<ConsensusMap>& input_maps, ConsensusMap& result_map)
  {
    if (input_maps.size() != 1)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "exactly one input map required");
    if (result_map.getColumnHeaders().size() != 2)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "two file descriptions required");
    if (result_map.getColumnHeaders().begin()->second.filename != result_map.getColumnHeaders().rbegin()->second.filename)
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "the two file descriptions have to contain the same file name");
    checkIds_(input_maps);

    //look up the light and heavy index
    Size light_index = numeric_limits<Size>::max();
    Size heavy_index = numeric_limits<Size>::max();
    for (ConsensusMap::ColumnHeaders::const_iterator it = result_map.getColumnHeaders().begin();
         it != result_map.getColumnHeaders().end();
         ++it)
    {
      if (it->second.label == "heavy")
      {
        heavy_index = it->first;
      }
      else if (it->second.label == "light")
      {
        light_index = it->first;
      }
    }
    if (light_index == numeric_limits<Size>::max() || heavy_index == numeric_limits<Size>::max())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "the input maps have to be labeled 'light' and 'heavy'");
    }

    result_map.clear(false);

    // sort consensus features by RT (and MZ) to speed up searching afterwards
    typedef ConstRefVector<ConsensusMap> RefMap;
    RefMap model_ref(input_maps[0].begin(), input_maps[0].end());
    model_ref.sortByPosition();

    //calculate matches
    ConsensusMap matches;
    //settings
    double rt_pair_dist = param_.getValue("rt_pair_dist");
    double rt_dev_low = param_.getValue("rt_dev_low");
    double rt_dev_high = param_.getValue("rt_dev_high");
    double mz_dev = param_.getValue("mz_dev");
    DoubleList mz_pair_dists = param_.getValue("mz_pair_dists");
    bool mrm = param_.getValue("mrm").toBool();

    //estimate RT parameters
    if (param_.getValue("rt_estimate") == "true")
    {
      //find all possible RT distances of features with the same charge and a good m/z distance
      vector<double> dists;
      dists.reserve(model_ref.size());
      for (RefMap::const_iterator it = model_ref.begin(); it != model_ref.end(); ++it)
      {
        for (RefMap::const_iterator it2 = model_ref.begin(); it2 != model_ref.end(); ++it2)
        {
          for (DoubleList::const_iterator dist_it = mz_pair_dists.begin(); dist_it != mz_pair_dists.end(); ++dist_it)
          {
            double mz_pair_dist = *dist_it;
            if (it2->getCharge() == it->getCharge()
               && it2->getMZ() >= it->getMZ() + mz_pair_dist / it->getCharge() - mz_dev
               && it2->getMZ() <= it->getMZ() + mz_pair_dist / it->getCharge() + mz_dev)
            {
              dists.push_back(it2->getRT() - it->getRT());
            }
          }
        }
      }
      if (dists.empty())
      {
        cout << "Warning: Could not find pairs for RT distance estimation. The manual settings are used!" << endl;
      }
      else
      {
        if (dists.size() < 50)
        {
          cout << "Warning: Found only " << dists.size() << " pairs. The estimated shift and std deviation are probably not reliable!" << endl;
        }
        //--------------------------- estimate initial parameters of fit ---------------------------
        GaussFitter::GaussFitResult result(-1, -1, -1);
        //first estimate of the optimal shift: median of the distances
        sort(dists.begin(), dists.end());
        Size median_index = dists.size() / 2;
        result.x0 = dists[median_index];
        //create histogram of distances
        //consider only the maximum of pairs, centered around the optimal shift
        Size max_pairs = model_ref.size() / 2;
        Size start_index = (Size) max((SignedSize)0, (SignedSize)(median_index - max_pairs / 2));
        Size end_index = (Size) min((SignedSize)(dists.size() - 1), (SignedSize)(median_index + max_pairs / 2));
        double start_value = dists[start_index];
        double end_value = dists[end_index];
        double bin_step = fabs(end_value - start_value) / 99.999; //ensure that we have 100 bins
        Math::Histogram<> hist(start_value, end_value, bin_step);
        //std::cout << "HIST from " << start_value << " to " << end_value << " (bin size " << bin_step << ")" << endl;
        for (Size i = start_index; i <= end_index; ++i)
        {
          hist.inc(dists[i]);
        }
        //cout << hist << endl;
        dists.clear();
        //determine median of bins (uniform background distribution)
        vector<Size> bins(hist.begin(), hist.end());
        sort(bins.begin(), bins.end());
        Size bin_median = bins[bins.size() / 2];
        bins.clear();
        //estimate scale A: maximum of the histogram
        Size max_value = hist.maxValue();
        result.A = max_value - bin_median;
        //overwrite estimate of x0 with the position of the highest bin
        for (Size i = 0; i < hist.size(); ++i)
        {
          if (hist[i] == max_value)
          {
            result.x0 = hist.centerOfBin(i);
            break;
          }
        }
        //estimate sigma: first time the count is less or equal the median count in the histogram
        double pos = result.x0;
        while (pos > start_value && hist.binValue(pos) > bin_median)
        {
          pos -= bin_step;
        }
        double sigma_low =  result.x0 - pos;
        pos = result.x0;
        while (pos<end_value&& hist.binValue(pos)> bin_median)
        {
          pos += bin_step;
        }
        double sigma_high = pos - result.x0;
        result.sigma = (sigma_high + sigma_low) / 6.0;
        //cout << "estimated optimal RT distance (before fit): " << result.x0 << endl;
        //cout << "estimated allowed deviation (before fit): " << result.sigma*3.0 << endl;
        //--------------------------- do gauss fit ---------------------------
        vector<DPosition<2> > points(hist.size());
        for (Size i = 0; i < hist.size(); ++i)
        {
          points[i][0] = hist.centerOfBin(i);
          points[i][1] = hist[i];
        }
        GaussFitter fitter;
        fitter.setInitialParameters(result);
        result = fitter.fit(points);
        cout << "estimated optimal RT distance: " << result.x0 << endl;
        cout << "estimated allowed deviation: " << fabs(result.sigma) * 3.0 << endl;
        rt_pair_dist = result.x0;
        rt_dev_low = fabs(result.sigma) * 3.0;
        rt_dev_high = fabs(result.sigma) * 3.0;
      }
    }


    // check each feature
    for (RefMap::const_iterator it = model_ref.begin(); it != model_ref.end(); ++it)
    {
      for (DoubleList::const_iterator dist_it = mz_pair_dists.begin(); dist_it != mz_pair_dists.end(); ++dist_it)
      {
        double mz_pair_dist = *dist_it;
        RefMap::const_iterator it2 = lower_bound(model_ref.begin(), model_ref.end(), it->getRT() + rt_pair_dist - rt_dev_low, ConsensusFeature::RTLess());
        while (it2 != model_ref.end() && it2->getRT() <= it->getRT() + rt_pair_dist + rt_dev_high)
        {
          // if in mrm mode, we need to compare precursor mass difference and fragment mass difference, charge remains the same

          double prec_mz_diff(0);
          if (mrm)
          {
            prec_mz_diff = fabs((double)it2->getMetaValue("MZ") - (double)it->getMetaValue("MZ"));
            if (it->getCharge() != 0)
            {
              prec_mz_diff = fabs(prec_mz_diff - mz_pair_dist / it->getCharge());
            }
            else
            {
              prec_mz_diff = fabs(prec_mz_diff - mz_pair_dist);
            }
          }

          bool mrm_correct_dist(false);
          double frag_mz_diff = fabs(it->getMZ() - it2->getMZ());

          //cerr << it->getRT() << " charge1=" << it->getCharge() << ", charge2=" << it2->getCharge() << ", prec_diff=" << prec_mz_diff << ", frag_diff=" << frag_mz_diff << endl;

          if (mrm &&
              it2->getCharge() == it->getCharge() &&
              prec_mz_diff < mz_dev &&
              (frag_mz_diff < mz_dev || fabs(frag_mz_diff - mz_pair_dist) < mz_dev))
          {
            mrm_correct_dist = true;
            //cerr << "mrm_correct_dist" << endl;
          }

          if ((mrm && mrm_correct_dist) || (!mrm &&
                                            it2->getCharge() == it->getCharge() &&
                                            it2->getMZ() >= it->getMZ() + mz_pair_dist / it->getCharge() - mz_dev &&
                                            it2->getMZ() <= it->getMZ() + mz_pair_dist / it->getCharge() + mz_dev
                                            ))
          {
            //cerr << "dist correct" << endl;
            double score = sqrt(
              PValue_(it2->getMZ() - it->getMZ(), mz_pair_dist / it->getCharge(), mz_dev, mz_dev) *
              PValue_(it2->getRT() - it->getRT(), rt_pair_dist, rt_dev_low, rt_dev_high)
              );

            // Note: we used to copy the id from the light feature here, but that strategy does not generalize to more than two labels.
            // We might want to report consensus features where the light one is missing but more than one heavier variant was found.
            // Also, the old strategy is inconsistent with what was done in the unlabeled case.  Thus now we assign a new unique id here.
            matches.push_back(ConsensusFeature());
            matches.back().setUniqueId();

            matches.back().insert(light_index, *it);
            matches.back().clearMetaInfo();
            matches.back().insert(heavy_index, *it2);
            matches.back().setQuality(score);
            matches.back().setCharge(it->getCharge());
            matches.back().computeMonoisotopicConsensus();
          }
          ++it2;
        }
      }
    }

    //compute best pairs
    // - sort matches by quality
    // - take highest-quality matches first (greedy) and mark them as used
    set<Size> used_features;
    matches.sortByQuality(true);
    for (ConsensusMap::const_iterator match = matches.begin(); match != matches.end(); ++match)
    {
      //check if features are not used yet
      if (used_features.find(match->begin()->getUniqueId()) == used_features.end() &&
          used_features.find(match->rbegin()->getUniqueId()) == used_features.end()
          )
      {
        //if unused, add it to the final set of elements
        result_map.push_back(*match);
        used_features.insert(match->begin()->getUniqueId());
        used_features.insert(match->rbegin()->getUniqueId());
      }
    }

    //Add protein identifications to result map
    for (Size i = 0; i < input_maps.size(); ++i)
    {
      result_map.getProteinIdentifications().insert(result_map.getProteinIdentifications().end(), input_maps[i].getProteinIdentifications().begin(), input_maps[i].getProteinIdentifications().end());
    }

    //Add unassigned peptide identifications to result map
    for (Size i = 0; i < input_maps.size(); ++i)
    {
      result_map.getUnassignedPeptideIdentifications().insert(result_map.getUnassignedPeptideIdentifications().end(), input_maps[i].getUnassignedPeptideIdentifications().begin(), input_maps[i].getUnassignedPeptideIdentifications().end());
    }

    // Very useful for checking the results, and the ids have no real meaning anyway
    result_map.sortByMZ();
  }

}
