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
		// defaults_.setValue("symmetric_regression","true","If true, linear regression will be based on (y-x) versus (x+y).\nIf false, a \"standard\" linear regression will be performed for y versus x.");
		// defaults_.setValidStrings("symmetric_regression",StringList::create("true,false"));
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
		Size reference_index = reference_index_ - 1; // local index is 0-based
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

/*
		startProgress(0, 10 * maps.size(), "aligning feature maps");

    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
    vector<ConsensusMap> input(2);
		ConsensusMap::convert(reference_index, maps[reference_index], input[0]);

		// init superimposer and pairfinder with model and parameters
		PoseClusteringAffineSuperimposer superimposer;
    superimposer.setParameters(param_.copy("superimposer:",true));
		superimposer.setLogType(getLogType());

    StablePairFinder pairfinder;
    pairfinder.setParameters(param_.copy("pairfinder:",true));
		pairfinder.setLogType(getLogType());

    for (Size i = 0; i < maps.size(); ++i)
		{
			setProgress(10*i);
			if (i != reference_index)
			{
				ConsensusMap::convert(i,maps[i],input[1]);
				setProgress(10*i+1);

				// run superimposer to find the global transformation
	      vector<TransformationDescription> si_trafos;
	      superimposer.run(input, si_trafos);
				setProgress(10*i+2);

				// apply transformation
				for (Size j=0; j<input[1].size(); ++j)
				{
					DoubleReal rt = input[1][j].getRT();
					si_trafos[0].apply(rt);
					input[1][j].setRT(rt);
					input[1][j].begin()->asMutable().setRT(rt);
				}
				setProgress(10*i+3);

	      //run pairfinder to find pairs
				ConsensusMap result;
				pairfinder.run(input, result);
				setProgress(10*i+4);

				// calculate the small local transformation
				TransformationDescription trafo;
				try
				{
					trafo = calculateRegression_(i,reference_index,result,symmetric_regression);
				}
				catch (Exception::Precondition& exception) // TODO is there a better way to deal with this situation?
				{
					cerr << "Warning: MapAlignementAlgorithmPoseClustering could not compute a refined mapping. Using initial estimation." << endl;
					trafo.setName("linear");
					trafo.setParam("slope",1.0);
					trafo.setParam("intercept",0.0);
				}
				setProgress(10*i+5);

				// combine the two transformations
				transformations[i].setName("linear");
				transformations[i].setParam("slope",(DoubleReal)si_trafos[0].getParam("slope")*(DoubleReal)trafo.getParam("slope"));
				transformations[i].setParam("intercept",(DoubleReal)trafo.getParam("slope")*(DoubleReal)si_trafos[0].getParam("intercept")+(DoubleReal)trafo.getParam("intercept"));

				// apply transformation (global and local)
				transformSingleFeatureMap(maps[i],transformations[i]);

				setProgress(10*i+6);
			}
		}
		if (reference_file_.empty())
		{
			//set no transformation for reference map:
			transformations[reference_index].setName("none");
		}
		else maps.resize(maps.size() - 1);

		setProgress(10*maps.size());
		endProgress();
	}
*/

/*
	TransformationDescription MapAlignmentAlgorithmPoseClustering::calculateRegression_(Size const index_x_map, Size const index_y_map, ConsensusMap const& consensus_map, bool const symmetric_regression) const
	{
		// the result
		TransformationDescription lm;

		// coordinates of appropriate pairs will be stored here
		vector<double> vec_x;
		vector<double> vec_y;

		// search for pairs, optionally apply coordinate transformation, store coordinates in vec_x and vec_y
		FeatureHandle probe_x(index_x_map,Peak2D(),0);
		FeatureHandle probe_y(index_y_map,Peak2D(),0);
		for ( ConsensusMap::const_iterator iter = consensus_map.begin(); iter != consensus_map.end(); ++iter )
		{
			DoubleReal rt_x, rt_y;
			ConsensusFeature::HandleSetType::const_iterator handle_x = iter->lower_bound(probe_x);
			if ( handle_x == iter->end() || handle_x->getMapIndex() != index_x_map)
			{
				continue;
			}
			rt_x = handle_x->getRT();

			ConsensusFeature::HandleSetType::const_iterator handle_y = iter->lower_bound(probe_y);
			if ( handle_y == iter->end() || handle_y->getMapIndex() != index_y_map)
			{
				continue;
			}
			rt_y = handle_y->getRT();

			if (symmetric_regression)
			{
				vec_x.push_back(rt_y+rt_x);
				vec_y.push_back(rt_y-rt_x);
			}
			else
			{
				vec_x.push_back(rt_x);
				vec_y.push_back(rt_y);
			}
		}

		// calculate the linear regression
		double slope, intercept;
		if ( vec_x.size() >= 2 ) // you would expect that
		{
			double cov00, cov01, cov11, sumsq;
			gsl_fit_linear(&vec_x.front(), 1, &vec_y.front(), 1, vec_x.size(), &intercept, &slope, &cov00, &cov01, &cov11, &sumsq);
		}
		else
		{
		}

		// assign final result, undo optional coordinate transform
		if (symmetric_regression)
		{
			lm.setName("linear");
			lm.setParam("slope", (1.+slope)/(1.-slope) );
			lm.setParam("intercept", intercept*1.41421356237309504880 ); // intercept*sqrt(2)
		}
		else
		{
			lm.setName("linear");
			lm.setParam("slope",slope);
			lm.setParam("intercept",intercept);
		}

		return lm;
	}
*/

	void MapAlignmentAlgorithmPoseClustering::getDefaultModel(String& model_type, 
																														Param& params)
	{
		model_type = "linear";
		params.clear();
		params.setValue("symmetric_regression", "true");
	}


} //namespace
