// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

	MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering()
		: MapAlignmentAlgorithm(), reference_index_(0), reference_file_()
	{
		setName("MapAlignmentAlgorithmPoseClustering");

		defaults_.insert("superimposer:",PoseClusteringAffineSuperimposer().getParameters());
		defaults_.insert("pairfinder:",StablePairFinder().getParameters());
		defaults_.setValue("max_num_peaks_considered",400,"The maximal number of peaks to be considered per map.  This cutoff is only applied to peak maps. To use all peaks, set this to '-1'.");
		defaults_.setMinInt("max_num_peaks_considered", -1);
		//TODO 'max_num_peaks_considered' should apply to peaks and features!! (Clemens)
		defaultsToParam_();
	}

	MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering()
	{
	}

	void MapAlignmentAlgorithmPoseClustering::setReference(Size reference_index, const String& reference_file)
	{
		reference_index_ = reference_index;
		// can't load the file yet because we don't know if the type will match:
		reference_file_ = reference_file;
	}

	void MapAlignmentAlgorithmPoseClustering::alignPeakMaps(vector< MSExperiment<> >& maps, vector<TransformationDescription>& transformations)
	{
		// prepare transformations for output
		transformations.clear();

		const Int max_num_peaks_considered = param_.getValue("max_num_peaks_considered");

		// reference map:
		Size reference_index = reference_index_ - 1; // local index is 0-based
		if (!reference_file_.empty())
		{
			if (FileHandler::getType(reference_file_) != FileTypes::MZML)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "reference file must be of type mzML in this case (same as input)");
			}
			maps.resize(maps.size() + 1);
			MzMLFile().load(reference_file_, maps.back());
			reference_index = maps.size() - 1;
		}
		else if (reference_index_ == 0) // no reference given
		{
			// use map with highest number of peaks as reference:
			Size max_count = 0;
			for (Size m = 0; m < maps.size(); ++m)
			{
				// initialize "getSize()" by calling "updateRanges()"
				maps[m].updateRanges(1);
				if (maps[m].getSize() > max_count)
				{
					max_count = maps[m].getSize();
					reference_index = m;
				}
			}
		}

		computeTransformations_(maps, transformations, reference_index, 
														max_num_peaks_considered);
	}


	template <typename MapType>
	void MapAlignmentAlgorithmPoseClustering::computeTransformations_(
		vector<MapType>& maps, vector<TransformationDescription>& transformations,
		Size reference_index, Size max_num_peaks_considered)
	{
		startProgress(0, 10 * maps.size(), "aligning input maps");

    // build a consensus map of the elements of the reference map (for consensus
    // map input, take only the highest peaks)
    vector<ConsensusMap> input(2);
		ConsensusMap::convert(reference_index, maps[reference_index], input[0], 
													max_num_peaks_considered);

		//init superimposer and pairfinder with model and parameters
		PoseClusteringAffineSuperimposer superimposer;
    superimposer.setParameters(param_.copy("superimposer:", true));
		superimposer.setLogType(getLogType());

    StablePairFinder pairfinder;
    pairfinder.setParameters(param_.copy("pairfinder:", true));
		pairfinder.setLogType(getLogType());

		for (Size i = 0; i < maps.size(); ++i)
		{
			setProgress(10 * i);
			if (i != reference_index)
			{
				// build scene_map
				ConsensusMap::convert(i, maps[i], input[1], max_num_peaks_considered);
				setProgress(10 * i + 1);

				// run superimposer to find the global transformation
	      vector<TransformationDescription> si_trafos;
	      superimposer.run(input, si_trafos);
				setProgress(10 * i + 2);

				// apply transformation to consensus features and contained feature
				// handles
				for (Size j = 0; j < input[1].size(); ++j)
				{
					//Calculate new RT
					DoubleReal rt = input[1][j].getRT();
					rt = si_trafos[0].apply(rt);
					//Set RT of consensus feature centroid
					input[1][j].setRT(rt);
					//Set RT of consensus feature handles
					input[1][j].begin()->asMutable().setRT(rt);
				}
				setProgress(10 * i + 3);

	      //run pairfinder fo find pairs
				ConsensusMap result;
				pairfinder.run(input, result);
				setProgress(10 * i + 4);

				// calculate the local transformation
				si_trafos[0].invert(); // to undo the transformation applied above
				TransformationDescription::DataPoints data;
				for (ConsensusMap::Iterator it = result.begin(); it != result.end();
						 ++it)
				{
					if (it->size() == 2) // two matching features
					{
						DoubleReal x = 0, y = 0;
						for (ConsensusFeature::iterator feat_it = it->begin();
								 feat_it != it->end(); ++feat_it)
						{
							// one feature should be from the reference map:
							if (feat_it->getMapIndex() == reference_index)
							{
								y = feat_it->getRT();
							}
							// transform RT back to the original scale:
							else x = si_trafos[0].apply(feat_it->getRT());
						}
						data.push_back(make_pair(x, y));
					}
				}
				setProgress(10 * i + 5);
				TransformationDescription trafo(data);
				transformations.push_back(trafo);
				setProgress(10*i+6);
			}

			else if (reference_file_.empty())
			{
				// set no transformation for reference map:
				TransformationDescription trafo;
				trafo.fitModel("identity");
				transformations.push_back(trafo);
			}
		}

		setProgress(10 * maps.size());
		endProgress();

		// reference file was added to "maps", has to be removed now:
		if (!reference_file_.empty()) maps.resize(maps.size() - 1);
	}


	void MapAlignmentAlgorithmPoseClustering::alignFeatureMaps(vector< FeatureMap<> >& maps, vector<TransformationDescription>& transformations)
	{
		// prepare transformations for output
		transformations.clear();

		// reference map:
		Size reference_index(0);
		if (!reference_file_.empty())
		{
			if (FileHandler::getType(reference_file_) != FileTypes::FEATUREXML)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "reference file must be of type featureXML in this case (same as input)");
			}
			maps.resize(maps.size() + 1);
			FeatureXMLFile().load(reference_file_, maps.back());
			reference_index = maps.size() - 1;
		}
    else if (reference_index_ == 0) // no reference given
    {
			// use map with highest number of features as reference:
      Size max_count = 0;
      for (Size m = 0; m < maps.size(); ++m)
      {
        if (maps[m].size() > max_count)
        {
          max_count = maps[m].size();
          reference_index = m;
        }
      }
    }

		computeTransformations_(maps, transformations, reference_index);
	}


	void MapAlignmentAlgorithmPoseClustering::getDefaultModel(String& model_type, 
																														Param& params)
	{
		model_type = "linear";
		params.clear();
		params.setValue("symmetric_regression", "true");
	}


} //namespace
