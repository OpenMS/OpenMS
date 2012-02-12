// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <cmath> // for "abs"
#include <limits> // for "max"

using namespace std;

namespace OpenMS
{

	MapAlignmentAlgorithmIdentification::MapAlignmentAlgorithmIdentification()
		: MapAlignmentAlgorithm(), reference_index_(0), reference_(),
			score_threshold_(0.0), min_run_occur_(0)
	{
		setName("MapAlignmentAlgorithmIdentification");

		defaults_.setValue("peptide_score_threshold", 0.0, "Score threshold for peptide hits to be used in the alignment.\nSelect a value that allows only 'high confidence' matches.");

		defaults_.setValue("min_run_occur", 2, "Minimum number of runs (incl. reference, if any) a peptide must occur in to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.");
		defaults_.setMinInt("min_run_occur", 2);

		defaults_.setValue("max_rt_shift", 0.5, "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale.");
		defaults_.setMinFloat("max_rt_shift", 0.0);
		
		defaults_.setValue("use_unassigned_peptides", "true", "Should unassigned peptide identifications be used when computing an alignment of feature maps? If 'false', only peptide IDs assigned to features will be used.");
		defaults_.setValidStrings("use_unassigned_peptides",
															StringList::create("true,false"));

		defaults_.setValue("use_feature_rt", "false", "When aligning feature maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'.");
		defaults_.setValidStrings("use_feature_rt", 
															StringList::create("true,false"));
		defaultsToParam_();
	}

	
	MapAlignmentAlgorithmIdentification::~MapAlignmentAlgorithmIdentification()
	{
	}


	void MapAlignmentAlgorithmIdentification::setReference(
		Size reference_index, const String& reference_file)
	{
		reference_.clear();
		reference_index_ = reference_index;
		// reference is one of the input files, or no reference given:
		if (reference_index_ || reference_file.empty()) return;

		// reference is external file:
		LOG_DEBUG << "Extracting reference RT data..." << endl;
		SeqToList rt_data;
		bool sorted = true;
		FileTypes::Type filetype = FileHandler::getType(reference_file);
		if (filetype == FileTypes::MZML)
		{
			MSExperiment<> experiment;
			MzMLFile().load(reference_file, experiment);
			getRetentionTimes_(experiment, rt_data);
			sorted = false;
		}
		else if (filetype == FileTypes::FEATUREXML)
		{
			FeatureMap<> features;
			FeatureXMLFile().load(reference_file, features);
			getRetentionTimes_(features, rt_data);
		}
		else if (filetype == FileTypes::CONSENSUSXML)
		{
			ConsensusMap features;
			ConsensusXMLFile().load(reference_file, features);
			getRetentionTimes_(features, rt_data);
		}
		else if (filetype == FileTypes::IDXML)
		{
			vector<ProteinIdentification> proteins;
			vector<PeptideIdentification> peptides;
			IdXMLFile().load(reference_file, proteins, peptides);
			getRetentionTimes_(peptides, rt_data);
		}

		computeMedians_(rt_data, reference_, sorted);
		if (reference_.empty())
		{
			throw Exception::MissingInformation(
				__FILE__, __LINE__, __PRETTY_FUNCTION__, 
				"Could not extract retention time information from the reference file");
		}
	}


	void MapAlignmentAlgorithmIdentification::checkParameters_(Size runs)
	{
		min_run_occur_ = param_.getValue("min_run_occur");
		// reference is not counted as a regular run:
		if (!reference_.empty()) runs++;
		if (min_run_occur_ > runs)
		{
			LOG_WARN << "Warning: Value of parameter 'min_run_occur' (here: " + String(min_run_occur_) + ") is higher than the number of runs incl. reference (here: " + String(runs) + "). Using " + String(runs) + " instead." << endl;
			min_run_occur_ = runs;
		}
		
		score_threshold_ = param_.getValue("peptide_score_threshold");		
	}


	void MapAlignmentAlgorithmIdentification::alignPeakMaps(
		vector<MSExperiment<> >& maps,
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning peak maps");

		if (reference_index_) // reference is one of the input files
		{
			SeqToList rt_data;
			getRetentionTimes_(maps[reference_index_ - 1], rt_data);
			computeMedians_(rt_data, reference_, false);			
		}

		// one set of RT data for each input map, except reference:
		vector<SeqToList> rt_data(maps.size() - bool(reference_index_));
		for (Size i = 0, j = 0; i < maps.size(); ++i)
		{
			if (i == reference_index_ - 1) continue; // skip reference map, if any
			getRetentionTimes_(maps[i], rt_data[j++]);
		}
		setProgress(1);

		computeTransformations_(rt_data, transformations);
		setProgress(2);

		// transformPeakMaps(maps, transformations);
		setProgress(3);
		endProgress();
	}


	template <typename MapType>
	void MapAlignmentAlgorithmIdentification::alignFeatureOrConsensusMaps_(
		vector<MapType>& maps, vector<TransformationDescription>& transformations)
	{
		if (reference_index_) // reference is one of the input files
		{
			SeqToList rt_data;
			getRetentionTimes_(maps[reference_index_ - 1], rt_data);
			computeMedians_(rt_data, reference_, true);
		}

		// one set of RT data for each input map, except reference:
		vector<SeqToList> rt_data(maps.size() - bool(reference_index_));
		for (Size i = 0, j = 0; i < maps.size(); ++i)
		{
			if (i == reference_index_ - 1) continue; // skip reference map, if any
			getRetentionTimes_(maps[i], rt_data[j++]);
		}
		setProgress(1);

		computeTransformations_(rt_data, transformations, true);
		setProgress(2);
	}


	void MapAlignmentAlgorithmIdentification::alignFeatureMaps(
		vector<FeatureMap<> >& maps,
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning feature maps");

		alignFeatureOrConsensusMaps_(maps, transformations);

		// transformFeatureMaps(maps, transformations);
		setProgress(3);
		endProgress();
	}


	void MapAlignmentAlgorithmIdentification::alignConsensusMaps(
		vector<ConsensusMap>& maps, 
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning consensus maps");

		alignFeatureOrConsensusMaps_(maps, transformations);

		// transformConsensusMaps(maps, transformations);
		setProgress(3);
		endProgress();
	}


	void MapAlignmentAlgorithmIdentification::alignPeptideIdentifications(
		vector<vector<PeptideIdentification> >& maps,
		vector<TransformationDescription>& transformations)
	{
		checkParameters_(maps.size());
		startProgress(0, 3, "aligning peptide identifications");

		if (reference_index_) // reference is one of the input files
		{
			SeqToList rt_data;
			getRetentionTimes_(maps[reference_index_ - 1], rt_data);
			computeMedians_(rt_data, reference_, true);			
		}

		// one set of RT data for each input map, except reference:
		vector<SeqToList> rt_data(maps.size() - bool(reference_index_));
		for (Size i = 0, j = 0; i < maps.size(); ++i)
		{
			if (i == reference_index_ - 1) continue; // skip reference map, if any
			getRetentionTimes_(maps[i], rt_data[j++]);
		}
		setProgress(1);

		computeTransformations_(rt_data, transformations, true);
		setProgress(2);

		// transformPeptideIdentifications(maps, transformations);
		setProgress(3);
		endProgress();
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
			DoubleReal median = Math::median(rt_it->second.begin(), 
																			 rt_it->second.end(), sorted);
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
		if (peptide.isHigherScoreBetter()) return score >= score_threshold_;
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


	template <typename MapType>
	void MapAlignmentAlgorithmIdentification::getRetentionTimes_(
		MapType& features, SeqToList& rt_data)
	{
		bool use_feature_rt = param_.getValue("use_feature_rt").toBool();
		for (typename MapType::Iterator feat_it = features.begin();
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
		
		// filter RT data (remove peptides that elute in several fractions):
		// TODO
		
		// compute RT medians:
		LOG_DEBUG << "Computing RT medians..." << endl;
		vector<SeqToValue> medians_per_run(size);
		for (Size i = 0; i < size; ++i)
		{
			computeMedians_(rt_data[i], medians_per_run[i], sorted);
		}
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

		// get reference retention time scale: either directly from reference file,
		// or compute consensus time scale
		bool reference_given = !reference_.empty(); // reference file given
		if (reference_given)
		{
			// remove peptides that don't occur in enough runs:
			LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
			SeqToValue temp;
			SeqToValue::iterator pos = temp.begin(); // to prevent segfault below
			for (SeqToValue::iterator ref_it = reference_.begin(); 
					 ref_it != reference_.end(); ++ref_it)
			{
				SeqToList::iterator med_it = medians_per_seq.find(ref_it->first);
				if ((med_it != medians_per_seq.end()) && 
						(med_it->second.size() + 1 >= min_run_occur_))
				{
					temp.insert(pos, *ref_it);
					pos = --temp.end(); // would cause segfault if "temp" was empty
				}
			}
			temp.swap(reference_);
		}
		else // compute overall RT median per sequence (median of medians per run)
		{	
			LOG_DEBUG << "Computing overall RT medians per sequence..." << endl;

			// remove peptides that don't occur in enough runs (at least two):
			LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
			SeqToList temp;
			SeqToList::iterator pos = temp.begin(); // to prevent segfault below
			for (SeqToList::iterator med_it = medians_per_seq.begin();
					 med_it != medians_per_seq.end(); ++med_it)
			{
				if (med_it->second.size() >= min_run_occur_)
				{
					temp.insert(pos, *med_it);
					pos = --temp.end(); // would cause segfault if "temp" was empty
				}
			}
			temp.swap(medians_per_seq);
			computeMedians_(medians_per_seq, reference_);
		}

		DoubleReal max_rt_shift = param_.getValue("max_rt_shift");
		if (max_rt_shift == 0)
		{
			max_rt_shift = numeric_limits<DoubleReal>::max();
		}
		else if (max_rt_shift <= 1)
		{ // compute max. allowed shift from overall retention time range:
			DoubleReal rt_range, rt_min = reference_.begin()->second, 
				rt_max = rt_min;
			for (SeqToValue::iterator it = ++reference_.begin(); 
					 it != reference_.end(); ++it)
			{
				rt_min = min(rt_min, it->second);
				rt_max = max(rt_max, it->second);
			}
			rt_range = rt_max - rt_min;
			max_rt_shift *= rt_range;
		}
		LOG_DEBUG << "Max. allowed RT shift (in seconds): " << max_rt_shift << endl;

		// generate RT transformations:
		LOG_DEBUG << "Generating RT transformations..." << endl;
		LOG_INFO << "\nAlignment based on:" << endl; // diagnostic output
		for (Size i = 0, offset = 0; i < size + 1; ++i)
		{
			if (i == reference_index_ - 1)
			{
				// if one of the input maps was used as reference, it has been skipped 
				// so far - now we have to consider it again:
				TransformationDescription trafo;
				trafo.fitModel("identity");
				transforms.push_back(trafo);
				LOG_INFO << "- 0 data points for sample " << i + 1 << " (reference)\n";
				offset = 1;
			}
			if (i >= size) break;
			// to be useful for the alignment, a peptide sequence has to occur in the
			// current run ("medians_per_run[i]"), but also in at least one other run
			// ("medians_overall"):
			TransformationDescription::DataPoints data;
			for (SeqToValue::iterator med_it = medians_per_run[i].begin();
					 med_it != medians_per_run[i].end(); ++med_it)
			{
				SeqToValue::const_iterator pos = reference_.find(med_it->first);
				if ((pos != reference_.end()) && 
						(fabs(med_it->second - pos->second) <= max_rt_shift))
				{ // found, and satisfies "max_rt_shift" condition!
					data.push_back(make_pair(med_it->second, pos->second));
				}
			}
			transforms.push_back(TransformationDescription(data));
			LOG_INFO << "- " << data.size() << " data points for sample " 
							 << i + offset + 1 << "\n";
		}
		LOG_INFO << endl;

		if (!reference_given) reference_.clear(); // delete temporary reference
	}

} //namespace
