// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmIdentification.h>

namespace OpenMS
{

  FeatureGroupingAlgorithmIdentification::FeatureGroupingAlgorithmIdentification() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmIdentification");
    defaults_.setValue("xcorr_threshold", 1.2, "Peptide identifications having a SEQUEST XCorr score smaller than this are discarded.");
    defaults_.setValue("rt_stdev_threshold", 100.0, "Maximum allowed standard deviation of retention times within a group");
    defaults_.setValue("mz_stdev_threshold", 1.0, "Maximum allowed standard deviation of mass-to-charge within a group");
    defaultsToParam_();
  }

  FeatureGroupingAlgorithmIdentification::~FeatureGroupingAlgorithmIdentification()
  {
  }

  struct PepHit
  {
    Size pep_map_nr;
    Size pep_feature_nr; // Careful!!!! (Don't change the maps after the feature_nr is stored!!)
    Size pep_ident_nr;
    Size pep_hit_nr;
    AASequence pep_sequence;
    DoubleReal pep_rt;
    DoubleReal pep_mz;
    DoubleReal pep_xcorr;
    String pep_id_algorithm;
  };

  struct SortPepHit
  {
    bool
    operator ()(const PepHit& a, const PepHit& b) const
    {
      if ( a.pep_sequence != b.pep_sequence )
      {
        return a.pep_sequence < b.pep_sequence;
      }
      else
      {
        return a.pep_xcorr > b.pep_xcorr;
      }
    }
  };

  struct SortPepHitbyMap
  {
    bool
    operator ()(const PepHit& a, const PepHit& b) const
    {
      return a.pep_map_nr < b.pep_map_nr;
    }
  };

  void
  FeatureGroupingAlgorithmIdentification::group(const std::vector<FeatureMap<> >& maps, ConsensusMap& out)
  {
    // check that the number of maps is ok
    if ( maps.size() < 2 )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "At least two maps must be given!");
    }

    // get the params
    DoubleReal xcorr_threshold = param_.getValue("xcorr_threshold");
    DoubleReal rt_stdev_threshold = param_.getValue("rt_stdev_threshold");
    DoubleReal mz_stdev_threshold = param_.getValue("mz_stdev_threshold");

    /* In the first step we scan through all peptide identifications.
     * We disregard unreliable peptide identifications having a SEQUEST XCorr score less than 1.2.
     * We check whether the RT and the m/z value of the precursor ion lies within the convex hull of a feature.
     * In this case we assign the peptide identification to the feature.
     * Each feature can be annotated with many peptide identifications originating from many MS/MS scans within the experiment.
     * The peptide identifications are already assigned to the features.
     * I need to filter out the ones with an XCorr > 1.2
     */

    /* In the second step we filter the peptide annotations with respect to the retention times of the features they are assigned to.
     * If a peptide identification is assigned to two features with very different RTs in one map, it is likely that one or both features are falsely annotated.
     * This observation is used to filter out dubious identifications which otherwise might give rise to incorrect consensus features in the ground truth.
     * For each peptide identification, we compute the mean μ and standard deviation σ of the RT positions of the features to which it is assigned.
     * If σ > 100 s, then the identification is considered dubious and removed from all features.
     * Moreover, the identification is removed from all features, if any, whose RT positions deviate by more than 2σ from μ.
     * These filters are applied for each experiment separately.
     */

    std::vector<FeatureMap<> > feature_maps = maps; // copy maps, so that they can be changed. // ???? größte map als erste!?

    for ( Size i = 0; i < feature_maps.size(); ++i ) // feature maps
    {
      // 	feature_maps[i].sortByRT();

      std::map<AASequence, Int> abs; // absolute number of peptides with this sequence
      std::map<AASequence, DoubleReal> num_mean; // sum of rt
      std::map<AASequence, DoubleReal> mean; // mean = num_mean / abs
      std::map<AASequence, DoubleReal> num_var; // denominator of the variance
      std::map<AASequence, DoubleReal> stand_dev; // standard deviation = sqrt(num_var / (abs - 1) )

      // features
      for ( Size j = 0; j < feature_maps[i].size(); ++j )
      {
        for ( Size k = 0; k < feature_maps[i][j].getPeptideIdentifications().size(); ++k ) // peptide identifications
        {
          std::vector<PeptideHit> peptide_hits;
          for ( Size l = 0; l < feature_maps[i][j].getPeptideIdentifications()[k].getHits().size(); ++l ) // peptide hits
          {
            PeptideHit peptide_hit = feature_maps[i][j].getPeptideIdentifications()[k].getHits()[l];

            // first step: identifications having a XCorr smaller than the threshold are discarded
            if ( DoubleReal(peptide_hit.getMetaValue("XCorr")) > xcorr_threshold )
            {
              peptide_hit.setMetaValue("IDAlgorithm", (String) "true");
              peptide_hits.push_back(peptide_hit);

              // while stepping through the peptide identifications, the maps needed to calculate mean and standard deviation are filled
              ++abs[peptide_hit.getSequence()];
              num_mean[peptide_hit.getSequence()] += maps[i][j].getRT();
            }
            else
            {
              peptide_hit.setMetaValue("IDAlgorithm", (String) "discarded: XCorr too small");
              peptide_hits.push_back(peptide_hit);
            }
          }
          feature_maps[i][j].getPeptideIdentifications()[k].setHits(peptide_hits);
        }
      }

      // mean
      for ( std::map<AASequence, Int>::iterator iter = abs.begin(); iter != abs.end(); ++iter )
      {
        mean[iter->first] = num_mean[iter->first] / (DoubleReal) abs[iter->first];
      }

      // standard deviation: num_var
      for ( Size j = 0; j < feature_maps[i].size(); ++j ) // features
      {
        for ( Size k = 0; k < feature_maps[i][j].getPeptideIdentifications().size(); ++k ) // peptide identifications
        {
          for ( Size l = 0; l < feature_maps[i][j].getPeptideIdentifications()[k].getHits().size(); ++l ) // peptide hits
          {
            PeptideHit peptide_hit = feature_maps[i][j].getPeptideIdentifications()[k].getHits()[l];

            // if the hit is not yet discarded, compute the denominator of the variance
            if ( peptide_hit.getMetaValue("IDAlgorithm") == "true" )
            {
              if ( num_var.find(peptide_hit.getSequence()) != num_var.end() )
              {
                num_var[peptide_hit.getSequence()] += (feature_maps[i][j].getRT() - mean[peptide_hit.getSequence()]) * (feature_maps[i][j].getRT()
                    - mean[peptide_hit.getSequence()]);
              }
              else
              {
                num_var[peptide_hit.getSequence()] = (feature_maps[i][j].getRT() - mean[peptide_hit.getSequence()]) * (feature_maps[i][j].getRT()
                    - mean[peptide_hit.getSequence()]);
              }
            }
          }
        }
      }

      // standard deviation
      for ( std::map<AASequence, DoubleReal>::iterator iter = num_var.begin(); iter != num_var.end(); ++iter )
      {
        // if only one element exists, the stand. dev. is 0
        if ( abs[iter->first] == 1 )
        {
          stand_dev[iter->first] = 0.0;
        }
        else // otherwise standard_deviation = sqrt(var)
        {
          stand_dev[iter->first] = sqrt(num_var[iter->first] / (DoubleReal) (abs[iter->first] - 1));
        }
      }

      // features
      for ( Size j = 0; j < feature_maps[i].size(); ++j )
      {
        for ( Size k = 0; k < feature_maps[i][j].getPeptideIdentifications().size(); ++k ) // petide identifications
        {
          std::vector<PeptideHit> peptide_hits;
          for ( Size l = 0; l < feature_maps[i][j].getPeptideIdentifications()[k].getHits().size(); ++l ) // peptide hits
          {
            PeptideHit peptide_hit = feature_maps[i][j].getPeptideIdentifications()[k].getHits()[l];

            // if the hit is not yet discarded, test if standard deviation is to high or rt-mean_rt deviates by more than 2*standard deviation
            if ( peptide_hit.getMetaValue("IDAlgorithm") == "true" )
            {
              if ( stand_dev[peptide_hit.getSequence()] > rt_stdev_threshold )
              {
                peptide_hit.setMetaValue("IDAlgorithm", (String) "discarded: rt standard deviation too big");
                peptide_hits.push_back(peptide_hit);
              }
              else if ( (stand_dev[peptide_hit.getSequence()] != 0.0) && (fabs(feature_maps[i][j].getRT() - mean[peptide_hit.getSequence()]) > 2
                  * stand_dev[peptide_hit.getSequence()]) ) // MAGIC ALERT: factor 2 corresponds to about 95% quantile in a gaussian
              {
                peptide_hit.setMetaValue("IDAlgorithm", (String) "discarded: rt deviates by more than 2*standard deviation from mean");
                peptide_hits.push_back(peptide_hit);
              }
              else
              {
                peptide_hits.push_back(peptide_hit);
              }
            }
            else
            {
              peptide_hits.push_back(peptide_hit);
            }
          }
          feature_maps[i][j].getPeptideIdentifications()[k].setHits(peptide_hits);
        }
      }
    }

    /* In the third step we compute an initial list of consensus features, in which features with identical identifications are grouped across maps.
     * In the previous steps we have computed a set of associations between peptide identifications from MS/MS and LC-MS features.
     * The consensus features in our ground truth should have unique peptide identifications.
     * Therefore we start by compiling a complete list of all peptide identifications over all experiments.
     * Then we step through this list and for each identification we find the best-scoring features associated with it,
     * but at most one from each experiment, and add these features to the corresponding consensus feature.
     * In this way we maximize the sum of XCorr values for the peptide identifications in a consensus feature.
     * We discard dubious consensus feature whose m/z standard deviation is greater than 1.
     */

    // std::vector< PepHit > pep_hits;
    std::map<AASequence, std::vector<PepHit> > pep_hits_initial;

    // holds all the identifications I discard in the steps three to five, so that the reason can be assigned to the features as a meta value
    std::vector<PepHit> discarded_pephits;

    for ( Size i = 0; i < feature_maps.size(); ++i ) // feature maps
    {
      std::vector<PepHit> pep_hits_map;

      for ( Size j = 0; j < feature_maps[i].size(); ++j ) // features
      {
        for ( Size k = 0; k < feature_maps[i][j].getPeptideIdentifications().size(); ++k ) // petide identifications
        {
          for ( Size l = 0; l < feature_maps[i][j].getPeptideIdentifications()[k].getHits().size(); ++l ) // peptide hits
          {
            PeptideHit peptide_hit = feature_maps[i][j].getPeptideIdentifications()[k].getHits()[l];
            if ( peptide_hit.getMetaValue("IDAlgorithm") == "true" ) // put all hits that are not yet discarded in a vector
            {
              DoubleReal xcorr_pep;
              if ( peptide_hit.metaValueExists("XCorr") )
              {
                xcorr_pep = peptide_hit.getMetaValue("XCorr");
              }
              else
              {
                xcorr_pep = 0.0;
              }
              PepHit pep_hit =
                { i, j, k, l, peptide_hit.getSequence(), feature_maps[i][j].getRT(), feature_maps[i][j].getMZ(), xcorr_pep, "true" };
              pep_hits_map.push_back(pep_hit);
            }
          }
        }
      }

      std::sort(pep_hits_map.begin(), pep_hits_map.end(), SortPepHit()); // sort hits by sequence and XCorr

      // for every sequence only put the hit with the highest score in the initial list of peptide hits (i.e. consensus features)
      if ( !pep_hits_map.empty() )
      {
        pep_hits_initial[pep_hits_map[0].pep_sequence].push_back(pep_hits_map[0]);
        for ( std::vector<PepHit>::iterator iter = pep_hits_map.begin(); iter != pep_hits_map.end() - 1; ++iter )
        {
          if ( iter->pep_sequence != (iter + 1)->pep_sequence ) // only one feature from every map for every identification
          {
            pep_hits_initial[(iter + 1)->pep_sequence].push_back(*(iter + 1));
          }
          else
          {
            (iter + 1)->pep_id_algorithm = "discarded: this identification exists with a higher score in the same map.";
            discarded_pephits.push_back(*(iter + 1));
          }
        }
        pep_hits_map.clear();
      }
    }

    // compute m/z standard deviation and discard consensus features having one higher than the threshold
    for ( std::map<AASequence, std::vector<PepHit> >::iterator itermap = pep_hits_initial.begin(); itermap != pep_hits_initial.end(); )
    {
      Int abs = 0;
      DoubleReal num_mean_mz = 0.0;
      DoubleReal num_var_mz = 0.0;

      // in generateGoldStandard.py also discarded:
      // 	DoubleReal num_mean_rt = 0.0;
      // 	DoubleReal num_var_rt = 0.0;

      for ( std::vector<PepHit>::iterator itervec = itermap->second.begin(); itervec != itermap->second.end(); ++itervec )
      {
        ++abs;
        num_mean_mz += itervec->pep_mz;
        // 	num_mean_rt += itervec->pep_rt;
      }

      DoubleReal mean_mz = num_mean_mz / (DoubleReal) abs;
      // 	DoubleReal mean_rt = num_mean_rt / (DoubleReal)abs;

      for ( std::vector<PepHit>::iterator itervec = itermap->second.begin(); itervec != itermap->second.end(); ++itervec )
      {
        num_var_mz += (itervec->pep_mz - mean_mz) * (itervec->pep_mz - mean_mz);
        // 		num_var_rt += (itervec->pep_rt - mean_rt) * (itervec->pep_rt - mean_rt);
      }

      DoubleReal stand_dev_mz = sqrt(num_var_mz / (DoubleReal) (abs - 1));
      // 	DoubleReal stand_dev_rt = sqrt( num_var_rt / (DoubleReal)(abs - 1) );
      // 	bool dev_rt_too_big = false;

      /*	for(std::vector<PepHit>::iterator itervec = itermap->second.begin(); itervec != itermap->second.end(); ++itervec)
       {
       if(fabs(itervec->pep_rt - mean_rt) > 2 * stand_dev_rt)
       {
       dev_rt_too_big = true;
       }
       }	*/

      if ( stand_dev_mz > mz_stdev_threshold ) // discard the consensus features which have a higher m/z standard deviation
      {
        for ( std::vector<PepHit>::iterator itervec = itermap->second.begin(); itervec != itermap->second.end(); ++itervec )
        {
          itervec->pep_id_algorithm = "discarded: mz standard deviation too big";
          discarded_pephits.push_back(*itervec);
        }
        pep_hits_initial.erase(itermap++);
      }
      /*	else if(dev_rt_too_big)
       {
       for(std::vector<PepHit>::iterator itervec = itermap->second.begin(); itervec != itermap->second.end(); ++itervec)
       {
       itervec->pep_id_algorithm = "discarded: rt-standard deviation of c.f. was too big";
       discarded_pephits.push_back(*itervec);
       }
       pep_hits_initial.erase(itermap++);
       }*/
      else
      {
        ++itermap;
      }
    }

    /* Let the total XCorr score of a consensus feature be defined as the sum of XCorr values of all features contained in it.
     * After step three, it is possible that a feature is contained in different consensus features from the initial list.
     * In the fourth step we reduce the initial list such that each feature is contained in at most one consensus feature,
     * whose total score is the largest among all consensus features containing it.
     * We have developed a simple "greedy" strategy to achieve this goal.
     * The purified list of candidate consensus features is sorted in order of decreasing total score.
     * In each step we extract a consensus feature with maximum total XCorr score from the list.
     * This consensus feature is added to the consensus map, and all consensus features having a non-empty intersection with it are also removed from the list.
     * The process is iterated until no more consensus features can be found, i. e., the list has become empty. */

    std::multimap<DoubleReal, std::vector<PepHit> > pep_hits_xcorr; // consensus features sorted by total XCorr score
    std::vector<std::vector<PepHit> > pep_hits_max_xcorr; // list to add the c.f. with highest score

    for ( std::multimap<AASequence, std::vector<PepHit> >::iterator itermap = pep_hits_initial.begin(); itermap != pep_hits_initial.end(); ++itermap )
    {
      DoubleReal xcorr_sum = 0.0;

      for ( std::vector<PepHit>::iterator itervec = itermap->second.begin(); itervec != itermap->second.end(); ++itervec )
      {
        xcorr_sum += itervec->pep_xcorr;
      }

      // insert c.f. and the total score as key in multimap, so that the map is automatically sorted by those
      pep_hits_xcorr.insert(std::pair<DoubleReal, std::vector<PepHit> >(xcorr_sum, (itermap->second)));
    }

    while ( !pep_hits_xcorr.empty() )
    {
      if ( ((--pep_hits_xcorr.end())->second).size() > 1 )
      {
        for ( std::multimap<DoubleReal, std::vector<PepHit> >::iterator itermap = pep_hits_xcorr.begin(); itermap != --pep_hits_xcorr.end(); )
        {
          // take the c.f. with the biggest total score and remove all c.f having a not-empty intersection with it
          bool same_feature = false;
          Size i = 0;

          while ( i < ((--pep_hits_xcorr.end())->second).size() ) // und immer dran denken: .end() zeigt hinter die map/den vector!!!
          {
            Size j = 0;
            while ( j < (itermap->second).size() )
            {
              if ( ((--pep_hits_xcorr.end())->second)[i].pep_feature_nr == (itermap->second)[j].pep_feature_nr
                  && ((--pep_hits_xcorr.end())->second)[i].pep_map_nr == (itermap->second)[j].pep_map_nr )
              {
                i += ((--pep_hits_xcorr.end())->second).size();
                j += (itermap->second).size();
                same_feature = true;
              }
              ++j;
            }
            ++i;
          }
          if ( same_feature )
          {
            for ( std::vector<PepHit>::iterator itervec = (itermap->second).begin(); itervec != (itermap->second).end(); ++itervec )
            {
              itervec->pep_id_algorithm = "discarded: this feature was already inserted in another consensus feature with a better score";
              discarded_pephits.push_back(*itervec);
            }
            pep_hits_xcorr.erase(itermap++);
          }
          else
          {
            ++itermap;
          }
        }

        // what happens if i only discard the concerned features from the list and not the whole c.f.??
        // nothing differnet of course. might make sense, when there are more maps?!
        // bad idea: increases the runtime a lot!!

        /*	std::multimap< DoubleReal, std::vector< PepHit > > pep_hits_xcorr_help;
         for(std::multimap<DoubleReal, std::vector< PepHit > >::iterator itermap = pep_hits_xcorr.begin(); itermap != --pep_hits_xcorr.end(); ++itermap)
         {
         std::cout << itermap->first << "\t";
         DoubleReal tot_xcorr = itermap->first;

         for(std::vector< PepHit >::iterator iterend = ((--pep_hits_xcorr.end())->second).begin(); iterend < ((--pep_hits_xcorr.end())->second).end(); ++iterend)
         {
         for(std::vector< PepHit >::iterator iterpep = itermap->second.begin(); iterpep < itermap->second.end(); )
         {
         if(iterend->pep_feature_nr == iterpep->pep_feature_nr && iterend->pep_map_nr == iterpep->pep_map_nr)
         {
         tot_xcorr -= iterpep->pep_xcorr;
         itermap->second.erase(iterpep++);
         }
         else
         {
         ++iterpep;
         }
         }
         }
         if(!(itermap->second).empty())
         {
         pep_hits_xcorr_help.insert( std::pair<DoubleReal, std::vector< PepHit > >(tot_xcorr, itermap->second) );
         }
         } */
        pep_hits_max_xcorr.push_back((--pep_hits_xcorr.end())->second);
        pep_hits_xcorr.erase(--pep_hits_xcorr.end());
        // 	pep_hits_xcorr = pep_hits_xcorr_help;
      }
      else
      {
        pep_hits_xcorr.erase(--pep_hits_xcorr.end());
      }
    }

    /* In the fifth step, we apply a final filter for outliers and dubious identifications by comparing retention times across maps.
     * We calculate the RT sample variance within all consensus features in the consensus map
     * and discard consensus features whose standard deviation is greater than 2 times the sample standard deviation.
     * Since this filter relies upon RT information and hence bears the risk of introducing bias into the ground truth,
     * we confirmed that the removed consensus features are indeed outliers by visual inspection.
     */

    std::vector<DoubleReal> rt_diffs;
    DoubleReal diff_sum = 0.0;
    Int diff_abs = 0;

    for ( std::vector<std::vector<PepHit> >::iterator iterout = pep_hits_max_xcorr.begin(); iterout != pep_hits_max_xcorr.end(); ++iterout )
    {
      // mean and standard deviation for every c.f.
      // 	Int abs = 0; // denominator mean (and stand. dev.)
      // 	DoubleReal num_mean = 0.0; // numerator mean
      for ( std::vector<PepHit>::iterator iterin = (*iterout).begin(); iterin != (*iterout).end(); ++iterin )
      {
        for ( std::vector<PepHit>::iterator iterin_2 = (*iterout).begin(); iterin_2 != (*iterout).end(); ++iterin_2 )
        {
          if ( iterin != iterin_2 )
          {
            DoubleReal rt_diff = fabs(iterin->pep_rt - iterin_2->pep_rt);
            rt_diffs.push_back(rt_diff);
            diff_sum += rt_diff;
            ++diff_abs;
          }
        }
      }
    }

    DoubleReal mean_diff = diff_sum / (DoubleReal) diff_abs;
    DoubleReal rt_diff_std_num = 0.0;

    for ( std::vector<DoubleReal>::iterator iterdiff = rt_diffs.begin(); iterdiff != rt_diffs.end(); ++iterdiff )
    {
      rt_diff_std_num += (*iterdiff - mean_diff) * (*iterdiff - mean_diff);
    }

    DoubleReal rt_diff_std = sqrt(rt_diff_std_num / (DoubleReal) (--diff_abs));

    std::vector<std::vector<PepHit> > pep_hits_max;

    for ( std::vector<std::vector<PepHit> >::iterator iterout = pep_hits_max_xcorr.begin(); iterout != pep_hits_max_xcorr.end(); ++iterout )
    {
      // mean and standard deviation for every c.f.
      Int abs = 0; // denominator mean (and stand. dev.)
      DoubleReal num_mean = 0.0; // numerator mean
      for ( std::vector<PepHit>::iterator iterin = (*iterout).begin(); iterin != (*iterout).end(); ++iterin )
      {
        ++abs;
        num_mean += iterin->pep_rt;
      }
      DoubleReal mean = num_mean / (DoubleReal) abs;
      DoubleReal num_std = 0.0;
      for ( std::vector<PepHit>::iterator iterin = (*iterout).begin(); iterin != (*iterout).end(); ++iterin )
      {
        num_std += (iterin->pep_rt - mean) * (iterin->pep_rt - mean);
      }
      DoubleReal std = sqrt(num_std / (DoubleReal) (--abs));
      if ( std > 2 * rt_diff_std )
      {
        for ( std::vector<PepHit>::iterator itervec = (*iterout).begin(); itervec != (*iterout).end(); ++itervec )
        {
          itervec->pep_id_algorithm = "discarded: consensus feature has stdev greater than 2*sample_std_dev";
          discarded_pephits.push_back(*itervec);
        }
        // 	pep_hits_max_xcorr.erase(iterout++);
      }
      else
      {
        // 	++iterout;
        pep_hits_max.push_back(*iterout);
      }
    }

    /* sortier die aussortierten mit gründen wieder ein*/

    for ( std::vector<PepHit>::iterator iter = discarded_pephits.begin(); iter != discarded_pephits.end(); ++iter )
    {
      std::vector<PeptideHit> peptide_hits;
      for ( Size l = 0; l < feature_maps[iter->pep_map_nr][iter->pep_feature_nr].getPeptideIdentifications()[iter->pep_ident_nr].getHits().size(); ++l ) // peptide hits
      {
        PeptideHit peptide_hit = feature_maps[iter->pep_map_nr][iter->pep_feature_nr].getPeptideIdentifications()[iter->pep_ident_nr].getHits()[l];

        if ( l == iter->pep_hit_nr ) // ist das der hit, den ich im vektor habe?
        {
          peptide_hit.setMetaValue("IDAlgorithm", iter->pep_id_algorithm);
          peptide_hits.push_back(peptide_hit);

        }
        else
        {
          peptide_hits.push_back(peptide_hit);
        }
      }
      feature_maps[iter->pep_map_nr][iter->pep_feature_nr].getPeptideIdentifications()[iter->pep_ident_nr].setHits(peptide_hits);
    }

    /* Now build a consensus map and store this in consensusXML format. nö, 2. wird im hauptprogramm gemacht.*/
    ConsensusMap consensus_map_0;
    ConsensusMap::convert(0, feature_maps[0], consensus_map_0);

    for ( std::vector<std::vector<PepHit> >::iterator iterout = pep_hits_max.begin(); iterout != pep_hits_max.end(); ++iterout )
    {
      std::sort((*iterout).begin(), (*iterout).end(), SortPepHitbyMap());

      if ( (*iterout)[0].pep_map_nr == 0 )
      {
        for ( std::vector<PepHit>::iterator iterin = (*iterout).begin() + 1; iterin != (*iterout).end(); ++iterin )
        {
          // andere map, gleiche feature-nr, wie in erster (bzw. nullter) map, aber nicht beim insert machen, sondern in das richtige consensus feature füllen.
          // kann in featuremap ja woanders gewesen sein
          consensus_map_0[(*iterout)[0].pep_feature_nr].insert(iterin->pep_map_nr, iterin->pep_feature_nr,
              feature_maps[iterin->pep_map_nr][iterin->pep_feature_nr]);
          consensus_map_0[(*iterout)[0].pep_feature_nr].computeConsensus();
        }
        consensus_map_0[(*iterout)[0].pep_feature_nr].setMetaValue("IDAlgorithm_usedID", ((*iterout)[0].pep_sequence).toString());
      }
      else
      {
        ConsensusFeature consensus_feature;
        for ( std::vector<PepHit>::iterator iterin = (*iterout).begin(); iterin != (*iterout).end(); ++iterin )
        {
          // müsste gehen, weil fertiges consensus feature gepushbackt wird.
          consensus_feature.insert(iterin->pep_map_nr, iterin->pep_feature_nr, feature_maps[iterin->pep_map_nr][iterin->pep_feature_nr]);
          consensus_feature.computeConsensus();
        }
        consensus_feature.setMetaValue("IDAlgorithm_usedID", ((*iterout)[0].pep_sequence).toString());
        consensus_map_0.push_back(consensus_feature);
      }
    }

    for ( Size i = 1; i < feature_maps.size(); ++i ) // feature maps
    {
      // Add protein identifications to result map
      consensus_map_0.getProteinIdentifications().insert(consensus_map_0.getProteinIdentifications().end(),
          feature_maps[i].getProteinIdentifications().begin(), feature_maps[i].getProteinIdentifications().end());
      // hmm, benutzt dann nur den letzten identif.-run als referenz und die nummern der proteine aus denen der letzten hinzugefügten map.

      // Add unassigned peptide identifications to result map
      consensus_map_0.getUnassignedPeptideIdentifications().insert(consensus_map_0.getUnassignedPeptideIdentifications().end(),
          feature_maps[i].getUnassignedPeptideIdentifications().begin(), feature_maps[i].getUnassignedPeptideIdentifications().end());
    }

    // die singletons fehlen noch!!!

    // replace result with temporary map
    out.swap(consensus_map_0);
    // copy back the input maps (they have been deleted while swapping)
    // out.getFileDescriptions() = lala.getFileDescriptions();

    return;
  }

} // namespace OpenMS

