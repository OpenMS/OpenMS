// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <boost/regex.hpp>

using namespace std;

namespace OpenMS
{
  ConsensusMapNormalizerAlgorithmMedian::ConsensusMapNormalizerAlgorithmMedian() = default;

  ConsensusMapNormalizerAlgorithmMedian::~ConsensusMapNormalizerAlgorithmMedian() = default;

  Size ConsensusMapNormalizerAlgorithmMedian::computeMedians(const ConsensusMap & map, vector<double>& medians, const String& acc_filter, const String& desc_filter)
  {
    Size number_of_maps = map.getColumnHeaders().size();
    vector<vector<double> > feature_int(number_of_maps);
    medians.resize(number_of_maps);

    // get map with most features, reserve space for feature_int (unequal vector lengths, 0-features omitted)
    ConsensusMap::ColumnHeaders::const_iterator map_with_most_features = map.getColumnHeaders().find(0);
    UInt map_with_most_features_idx = 0;
    for (UInt i = 0; i < number_of_maps; i++)
    {
      ConsensusMap::ColumnHeaders::const_iterator it = map.getColumnHeaders().find(i);
      if (it == map.getColumnHeaders().end()) 
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(i));
      }
      else if (i >= feature_int.size())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String(i) + " exceeds map number");
      }
       
      feature_int[i].reserve(it->second.size);

      if (it->second.size > map_with_most_features->second.size)
      {
        map_with_most_features = it;
        map_with_most_features_idx = i;
      }
    }

    // fill feature_int with intensities
    Size pass_counter = 0;
    ConsensusMap::ConstIterator cf_it;
    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      if (!passesFilters_(cf_it, map, acc_filter, desc_filter))
      {
        continue;
      }
      ++pass_counter;

      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        feature_int[f_it->getMapIndex()].push_back(f_it->getIntensity());
      }
    }

    OPENMS_LOG_INFO << endl << "Using " << pass_counter << "/" << map.size() <<  " consensus features for computing normalization coefficients" << endl << endl;

    // do we have enough features passing the filters to compute the median for every map?
    bool enough_features_left = true;
    for (UInt j = 0; j < number_of_maps; j++)
    {
      //set all medians to 1.0 for now, so normalization will have no effect if we return
      medians[j] = 1.0;

      vector<double>& ints_j = feature_int[j];
      if (ints_j.empty())
      {
        enough_features_left = false;
      }
    }

    if (!enough_features_left)
    {
      OPENMS_LOG_WARN << endl << "Not enough features passing filters. Cannot compute normalization coefficients for all maps. Result will be unnormalized." << endl << endl;
      return 0;
    }
    else
    {
      //compute medians
      for (UInt j = 0; j < number_of_maps; j++)
      {
        vector<double>& ints_j = feature_int[j];
        medians[j] = Math::median(ints_j.begin(), ints_j.end());
      }
    }

    return map_with_most_features_idx;
  }

  void ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(ConsensusMap & map, NormalizationMethod method, const String& acc_filter, const String& desc_filter)
  {
    if (method == NM_SHIFT)
    {
      OPENMS_LOG_WARN << endl << "WARNING: normalization using median shifting is not recommended for regular log-normal MS data. Use this only if you know exactly what you're doing!" << endl << endl;
    }

    ConsensusMap::Iterator cf_it;
    ProgressLogger progresslogger;
    progresslogger.setLogType(ProgressLogger::CMD);
    progresslogger.startProgress(0, map.size(), "normalizing maps");

    vector<double> medians;
    Size index_of_largest_map = computeMedians(map, medians, acc_filter, desc_filter);

    for (cf_it = map.begin(); cf_it != map.end(); ++cf_it)
    {
      progresslogger.setProgress(cf_it - map.begin());
      ConsensusFeature::HandleSetType::const_iterator f_it;
      for (f_it = cf_it->getFeatures().begin(); f_it != cf_it->getFeatures().end(); ++f_it)
      {
        Size map_index = f_it->getMapIndex();
        if (method == NM_SCALE)
        {
          // scale to median of map with largest number of features
          f_it->asMutable().setIntensity(f_it->getIntensity() * medians[index_of_largest_map] / medians[map_index]);
        }
        else // method == NM_SHIFT
        {
          // shift to median of map with largest median in order to avoid negative intensities
          double max_median(numeric_limits<double>::min());
          Size max_median_index(0);
          for (Size i = 0; i < medians.size(); ++i)
          {
            if (medians[i] > max_median)
            {
              max_median = medians[i];
              max_median_index = i;
            }
          }
          f_it->asMutable().setIntensity(f_it->getIntensity() + medians[max_median_index] - medians[map_index]);
        }
      }
    }
    progresslogger.endProgress();
  }

  bool ConsensusMapNormalizerAlgorithmMedian::passesFilters_(ConsensusMap::ConstIterator cf_it, const ConsensusMap& map, const String& acc_filter, const String& desc_filter)
  {
    boost::regex acc_regexp(acc_filter);
    boost::regex desc_regexp(desc_filter);
    boost::cmatch m;

    if ((acc_filter.empty() || boost::regex_search("", m, acc_regexp)) &&
        (desc_filter.empty() || boost::regex_search("", m, desc_regexp)))
    {
      // feature passes (even if it has no identification!)
      return true;
    }

    const vector<ProteinIdentification>& prot_ids = map.getProteinIdentifications();
    const vector<PeptideIdentification>& pep_ids = cf_it->getPeptideIdentifications();

    for (vector<PeptideIdentification>::const_iterator p_it = pep_ids.begin(); p_it != pep_ids.end(); ++p_it)
    {
      const vector<PeptideHit>& hits = p_it->getHits();
      for (vector<PeptideHit>::const_iterator h_it = hits.begin(); h_it != hits.end(); ++h_it)
      {
        const set<String>& accs = h_it->extractProteinAccessionsSet();
        for (set<String>::const_iterator acc_it = accs.begin(); acc_it != accs.end(); ++acc_it)
        {
          // does accession match?
          if (!(acc_filter.empty() ||
                boost::regex_search("", m, acc_regexp) ||
                boost::regex_search(acc_it->c_str(), m, acc_regexp)))
          {
            //no
            continue;
          }

          // yes. does description match, too?
          if (desc_filter.empty() || boost::regex_search("", m, desc_regexp))
          {
            return true;
          }
          for (vector<ProteinIdentification>::const_iterator pr_it = prot_ids.begin(); pr_it != prot_ids.end(); ++pr_it)
          {
            std::vector<ProteinHit>::const_iterator pr_hit = const_cast<ProteinIdentification&>(*pr_it).findHit(*acc_it);
            if (pr_hit != pr_it->getHits().end())
            {
              const char* desc = pr_hit->getDescription().c_str();
              if (boost::regex_search(desc, m, desc_regexp))
              {
                return true;
              }
            }
          }
        }
      }
    }

    return false;
  }

}
