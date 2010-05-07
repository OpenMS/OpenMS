// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

#include <cmath> // for "abs"
#include <limits> // for "max"

using namespace std;

namespace OpenMS
{

	MapAlignmentAlgorithmIdentification::MapAlignmentAlgorithmIdentification()
		: MapAlignmentAlgorithm()
	{
		setName("MapAlignmentAlgorithmIdentification");

		defaults_.setValue("num_breakpoints", 10, "Number of breakpoints of the cubic spline in the smoothing step.\nThe breakpoints are spaced uniformly on the retention time interval.\nMore breakpoints mean less smoothing.");
		defaults_.setMinInt("num_breakpoints", 2);

		defaults_.setValue("peptide_score_threshold", 0.0, "Score threshold for peptide hits to be used in the alignment.\nSelect a value that allows only 'high confidence' matches.");

		defaults_.setValue("min_run_occur", 2, "Minimum number of runs a peptide must occur in to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.");
		defaults_.setMinInt("min_run_occur", 2);

		defaults_.setValue("max_rt_shift", 0.5, "Maximum realistic RT difference for a peptide (median per run vs. median overall).\nPeptides with higher shifts are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the total retention time range.");
		defaults_.setMinFloat("max_rt_shift", 0.0);
		
		defaults_.setValue("use_unassigned_peptides", "true", "Should unassigned peptide identifications be used when computing an alignment of feature maps?\nIf 'false', only peptide IDs assigned to features will be used.");
		defaults_.setValidStrings("use_unassigned_peptides",
															StringList::create("true,false"));

		defaults_.setValue("use_feature_rt", "false", "When aligning feature maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to.\nIf different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'.");
		defaults_.setValidStrings("use_feature_rt", 
															StringList::create("true,false"));
		defaultsToParam_();
	}

	
	MapAlignmentAlgorithmIdentification::~MapAlignmentAlgorithmIdentification()
	{
	}


	void MapAlignmentAlgorithmIdentification::checkParameters_(const Size runs)
	{
		Size min_run_occur = param_.getValue("min_run_occur");
		if (min_run_occur > runs)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Value of parameter 'min_run_occur' (here: " + String(min_run_occur) + ") must not exceed the number of runs (here: " + String(runs) + ")");
		}
	}


	void MapAlignmentAlgorithmIdentification::alignPeakMaps(
		vector<MSExperiment<> >& maps,
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning peak maps");

		score_threshold_ = param_.getValue("peptide_score_threshold");
		
		vector<SeqToList> rt_data(maps.size());
		for (Size i = 0; i < maps.size(); ++i)
		{
			getRetentionTimes_(maps[i], rt_data[i]);
		}
		setProgress(1);

		computeTransformations_(rt_data, transformations);
		setProgress(2);

		transformPeakMaps(maps, transformations);
		setProgress(3);
		endProgress();
	}


	void MapAlignmentAlgorithmIdentification::alignFeatureMaps(
		vector<FeatureMap<> >& maps,
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning feature maps");

		score_threshold_ = param_.getValue("peptide_score_threshold");
		
		vector<SeqToList> rt_data(maps.size());
		for (Size i = 0; i < maps.size(); ++i)
		{
			getRetentionTimes_(maps[i], rt_data[i]);
		}
		setProgress(1);

		computeTransformations_(rt_data, transformations, true);
		setProgress(2);

		transformFeatureMaps(maps, transformations);
		setProgress(3);
		endProgress();
	}


	void MapAlignmentAlgorithmIdentification::alignPeptideIdentifications(
		vector<vector<PeptideIdentification> >& maps,
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning peptide identifications");

		score_threshold_ = param_.getValue("peptide_score_threshold");
		
		vector<SeqToList> rt_data(maps.size());
		for (Size i = 0; i < maps.size(); ++i)
		{
			getRetentionTimes_(maps[i], rt_data[i]);
		}
		setProgress(1);

		computeTransformations_(rt_data, transformations, true);
		setProgress(2);

		transformPeptideIdentifications(maps, transformations);
		setProgress(3);
		endProgress();
	}
	
	
	// "values" will be sorted (unless "sorted" is true)
	DoubleReal MapAlignmentAlgorithmIdentification::median_(DoubleList& values,
																													bool sorted)
	{
		Size size = values.size();
		if (!size)
		{
			throw(Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
																			 "input sequence is empty"));
		}		
		if (!sorted) sort(values.begin(), values.end());
		if (size % 2 == 0) // even size => average two middle values
		{
			size /= 2;
			return (values[size - 1] + values[size]) / 2.0;
		}
		return values[(size - 1) / 2];
	}


	// RT lists in "rt_data" will be sorted (unless "sorted" is true)
	void MapAlignmentAlgorithmIdentification::computeMedians_(SeqToList& rt_data,
																														SeqToValue& medians,
																														bool sorted)
	{
		medians.clear();
		SeqToValue::iterator pos = medians.begin(); // prevent segfault (see below)
		for (SeqToList::iterator rt_it = rt_data.begin();
				 rt_it != rt_data.end(); ++rt_it)
		{
			DoubleReal median = median_(rt_it->second, sorted);
			medians.insert(pos,	make_pair(rt_it->first, median));
			pos = --medians.end(); // would cause segfault if "medians" were empty
		}
	}


	// list of peptide hits in "peptide" will be sorted
	bool MapAlignmentAlgorithmIdentification::hasGoodHit_(PeptideIdentification&
																												peptide)
	{
		if (peptide.empty() || peptide.getHits().empty())
			return false;
		peptide.sort();
		DoubleReal score = peptide.getHits().begin()->getScore();
		if (peptide.isHigherScoreBetter())
			return score >= score_threshold_;
		return score <= score_threshold_;
	}


	// lists of peptide hits in "peptides" will be sorted
	void MapAlignmentAlgorithmIdentification::getRetentionTimes_(
		vector<PeptideIdentification>& peptides, SeqToList& rt_data)
	{
		for (vector<PeptideIdentification>::iterator pep_it =
					 peptides.begin(); pep_it != peptides.end(); ++pep_it)
		{
			if (hasGoodHit_(*pep_it))
			{
				rt_data[pep_it->getHits()[0].getSequence().toString()] <<
					pep_it->getMetaValue("RT");
			}
		}	
	}

	
	// lists of peptide hits in "maps" will be sorted
	void MapAlignmentAlgorithmIdentification::getRetentionTimes_(
		MSExperiment<>& experiment, SeqToList& rt_data)
	{
		for (MSExperiment<>::Iterator exp_it = experiment.begin();
				 exp_it != experiment.end(); ++exp_it)
		{
			getRetentionTimes_(exp_it->getPeptideIdentifications(), rt_data);
		}
		// duplicates should not be possible -> no need to remove them
	}
	

	// lists of peptide hits in "features" will be sorted
	void MapAlignmentAlgorithmIdentification::getRetentionTimes_(
		FeatureMap<>& features, SeqToList& rt_data)
	{
		bool use_feature_rt = param_.getValue("use_feature_rt").toBool();
		for (FeatureMap<>::Iterator feat_it = features.begin();
				 feat_it != features.end(); ++feat_it)
		{
			if (use_feature_rt)
			{
				// find the peptide ID closest in RT to the feature centroid:
				String sequence;
				DoubleReal rt_distance = numeric_limits<DoubleReal>::max();
				bool any_good_hit = false;
				for (vector<PeptideIdentification>::iterator pep_it =
							 feat_it->getPeptideIdentifications().begin(); pep_it != 
							 feat_it->getPeptideIdentifications().end(); ++pep_it)
				{
					if (hasGoodHit_(*pep_it))
					{
						any_good_hit = true;
						DoubleReal current_distance = 
							abs(double(pep_it->getMetaValue("RT")) - feat_it->getRT());
						if (current_distance < rt_distance)
						{
							sequence = pep_it->getHits()[0].getSequence().toString();
							rt_distance = current_distance;
						}				
					}
				}
				if (any_good_hit) rt_data[sequence] << feat_it->getRT();
			}
			else
			{
				getRetentionTimes_(feat_it->getPeptideIdentifications(), rt_data);
			}
		}

		if (!use_feature_rt && param_.getValue("use_unassigned_peptides").toBool())
		{
			getRetentionTimes_(features.getUnassignedPeptideIdentifications(),
												 rt_data);
		}

		// remove duplicates (can occur if a peptide ID was assigned to several
		// features due to overlap or annotation tolerance):
		for (SeqToList::iterator rt_it = rt_data.begin();
				 rt_it != rt_data.end(); ++rt_it)
		{
			DoubleList& rt_values = rt_it->second;
			sort(rt_values.begin(), rt_values.end());
			DoubleList::iterator it = unique(rt_values.begin(), rt_values.end());
			rt_values.resize(it - rt_values.begin());
		}
	}


	void MapAlignmentAlgorithmIdentification::computeTransformations_(
		vector<SeqToList>& rt_data, vector<TransformationDescription>& transforms,
		bool sorted)
	{
		Size size = rt_data.size();
		transforms.clear();
		transforms.resize(size);
		
		// filter RT data (remove peptides that elude in several fractions):
		// TODO
		
		// compute RT medians (RT lists are already sorted):
		// cout << "Computing RT medians..." << endl;
		vector<SeqToValue> medians_per_run(size);
		for (Size i = 0; i < size; ++i)
		{
			computeMedians_(rt_data[i], medians_per_run[i], sorted);
		}

		// compute overall RT median per sequence (median of medians per run):
		// cout << "Computing overall RT medians per sequence..." << endl;
		SeqToList medians_per_seq;
		for (vector<SeqToValue>::iterator run_it = medians_per_run.begin();
				 run_it != medians_per_run.end(); ++run_it)
		{
			for (SeqToValue::iterator med_it = run_it->begin();
					 med_it != run_it->end(); ++med_it)
			{
				medians_per_seq[med_it->first] << med_it->second;
			}
		}
		// remove peptides that don't occur in enough runs (at least two):
		// cout << "Removing peptides that occur in too few runs..." << endl;
		SeqToList temp;
		SeqToList::iterator pos = temp.begin(); // to prevent segfault (see below)
		Size min_run_occur = param_.getValue("min_run_occur");
		for (SeqToList::iterator med_it = medians_per_seq.begin();
				 med_it != medians_per_seq.end(); ++med_it)
		{
			if (med_it->second.size() >= min_run_occur)
			{
				temp.insert(pos, *med_it);
				pos = --temp.end(); // would cause segfault if "temp" were empty
			}
		}
		temp.swap(medians_per_seq);
		SeqToValue medians_overall;
		computeMedians_(medians_per_seq, medians_overall);

		DoubleReal max_rt_shift = param_.getValue("max_rt_shift");
		if (max_rt_shift == 0)
		{
			max_rt_shift = numeric_limits<DoubleReal>::max();
		}
		else if (max_rt_shift <= 1)
		{ // compute max. allowed shift from overall retention time range:
			DoubleReal rt_range, rt_min = medians_overall.begin()->second, 
				rt_max = rt_min;
			for (SeqToValue::iterator it = ++medians_overall.begin(); 
					 it != medians_overall.end(); ++it)
			{
				rt_min = min(rt_min, it->second);
				rt_max = max(rt_max, it->second);
			}
			rt_range = rt_max - rt_min;
			max_rt_shift *= rt_range;
		}
			
		// generate RT transformations:
		// cout << "Generating RT transformations..." << endl;
		Int num_breakpoints = param_.getValue("num_breakpoints");
		LOG_INFO << "\nAlignment based on:" << endl; // diagnostic output
		for (Size i = 0; i < size; ++i)
		{
			TransformationDescription::PairVector pairs;
			// to be useful for the alignment, a peptide sequence has to occur in the
			// current run ("medians_per_run[i]"), but also in at least one other run
			// ("medians_overall")
			for (SeqToValue::iterator med_it = medians_per_run[i].begin();
					 med_it != medians_per_run[i].end(); ++med_it)
			{
				SeqToValue::const_iterator pos = medians_overall.find(med_it->first);
				if ((pos != medians_overall.end()) && 
						(fabs(med_it->second - pos->second) <= max_rt_shift))
				{ // found, and satisfies "max_rt_shift" condition!
					pairs.push_back(make_pair(med_it->second, pos->second));
				}
			}
		
			// maybe TODO: check "pairs.size()" here and abort if not enough data?
			transforms[i].setName("b_spline");
			transforms[i].setParam("num_breakpoints", num_breakpoints);
			transforms[i].setPairs(pairs);
			LOG_INFO << "- " << pairs.size() << " data points for sample " << i + 1
							 << endl;
		}
		LOG_INFO << endl;
	}

} //namespace
