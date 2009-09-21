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
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>

#include <gsl/gsl_fit.h>

#include <iostream>


namespace OpenMS
{

	MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering()
		: MapAlignmentAlgorithm()
	{
		setName("MapAlignmentAlgorithmPoseClustering");

		defaults_.insert("superimposer:",PoseClusteringAffineSuperimposer().getParameters());
		defaults_.insert("pairfinder:",StablePairFinder().getParameters());
		defaults_.setValue("symmetric_regression","true","If true, linear regression will be based on (y-x) versus (x+y).\nIf false, a \"standard\" linear regression will be performed for y versus x.");
		defaults_.setValidStrings("symmetric_regression",StringList::create("true,false"));
		defaults_.setValue("max_num_peaks_considered",400,"The maximal number of peaks to be considered per map.  This cutoff is only applied to peak maps.  For using all peaks, set this to -1.");
		defaults_.setMinInt("max_num_peaks_considered",-1);
		//TODO 'max_num_peaks_considered' should apply to peaks and features!! (Clemens)
		defaultsToParam_();
	}

	MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering()
	{
	}

	void MapAlignmentAlgorithmPoseClustering::alignPeakMaps(std::vector< MSExperiment<> >& maps, std::vector<TransformationDescription>& transformations)
	{
		// prepare transformations for output
		transformations.clear();
		transformations.resize(maps.size());

		const Int max_num_peaks_considered = param_.getValue("max_num_peaks_considered");
		const bool symmetric_regression = param_.getValue("symmetric_regression").toBool();

		//define reference map (the one with most peaks)
		Size reference_map_index = 0;
		Size max_count = 0;
		for (Size m=0; m<maps.size(); ++m)
		{
			// initialize getSize() by calling updateRanges()
			maps[m].updateRanges(1);
			if (maps[m].getSize()>max_count)
			{
				max_count = maps[m].getSize();
				reference_map_index = m;
			}
		}

    // build a consensus map of the elements of the reference map (take the 400 highest peaks)
    std::vector<ConsensusMap> input(2);
		ConsensusMap::convert( reference_map_index, maps[reference_map_index], input[0], max_num_peaks_considered );

		//init superimposer and pairfinder with model and parameters
		PoseClusteringAffineSuperimposer superimposer;
    superimposer.setParameters(param_.copy("superimposer:",true));

    StablePairFinder pairfinder;
    pairfinder.setParameters(param_.copy("pairfinder:",true));

		for (Size i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				// build scene_map
				ConsensusMap::convert( i, maps[i], input[1], max_num_peaks_considered );

				// run superimposer to find the global transformation
	      std::vector<TransformationDescription> si_trafos;
	      superimposer.run(input, si_trafos);

				//apply transformation to consensus feature and contained feature handles
				for (Size j=0; j<input[1].size(); ++j)
				{
					//Calculate new RT
					DoubleReal rt = input[1][j].getRT();
					si_trafos[0].apply(rt);
					//Set rt of consensus feature centroid
					input[1][j].setRT(rt);
					//Set RT of consensus feature handles
					input[1][j].begin()->asMutable().setRT(rt);
				}

	      //run pairfinder fo find pairs
				ConsensusMap result;
				pairfinder.run(input,result);

				// calculate the local transformation
				TransformationDescription trafo;
				try
				{
					trafo = calculateRegression_(i,reference_map_index,result,symmetric_regression);
				}
				catch (Exception::Precondition & /*exception*/ ) // TODO is there a better way to deal with this situation?
				{
					std::cerr << "Warning: MapAlignementAlgorithmPoseClustering could not compute a refined mapping. Using initial estimation." << std::endl;
					trafo.setName("linear");
					trafo.setParam("slope",1.0);
					trafo.setParam("intercept",0.0);
				}


				// combine the two transformations
				transformations[i].setName("linear");
				transformations[i].setParam("slope",(DoubleReal)si_trafos[0].getParam("slope")*(DoubleReal)trafo.getParam("slope"));
				transformations[i].setParam("intercept",(DoubleReal)trafo.getParam("slope")*(DoubleReal)si_trafos[0].getParam("intercept")+(DoubleReal)trafo.getParam("intercept"));

				// apply transformation to all scans
				for (Size j=0; j< maps[i].size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					transformations[i].apply(rt);
					maps[i][j].setRT(rt);
				}
			}
		}
		// set no transformation for reference map
		transformations[reference_map_index].setName("none");
	}


	void MapAlignmentAlgorithmPoseClustering::alignFeatureMaps(std::vector< FeatureMap<> >& maps, std::vector<TransformationDescription>& transformations)
	{
		// prepare transformations for output
		transformations.clear();
		transformations.resize(maps.size());

		startProgress(0, 10 * maps.size(),"aligning feature maps");

		const bool symmetric_regression = param_.getValue("symmetric_regression").toBool();

		// define reference map (the one with most peaks)
		Size reference_map_index = 0;
		Size max_count = 0;
		for (Size m=0; m<maps.size(); ++m)
		{
			if (maps[m].size()>max_count)
			{
				max_count = maps[m].size();
				reference_map_index = m;
			}
		}

    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
    std::vector<ConsensusMap> input(2);
		ConsensusMap::convert(reference_map_index, maps[reference_map_index], input[0]);

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
			if (i != reference_map_index)
			{
				ConsensusMap::convert(i,maps[i],input[1]);
				setProgress(10*i+1);

				// run superimposer to find the global transformation
	      std::vector<TransformationDescription> si_trafos;
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
					trafo = calculateRegression_(i,reference_map_index,result,symmetric_regression);
				}
				catch (Exception::Precondition& /*exception*/ ) // TODO is there a better way to deal with this situation?
				{
					std::cerr << "Warning: MapAlignementAlgorithmPoseClustering could not compute a refined mapping. Using initial estimation." << std::endl;
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
#if 1 // do it the new way - "deep"
				transformSingleFeatureMap(maps[i],transformations[i]);
#else // do it the old way - "shallow", does not transform convex hulls etc.! - TODO remove this piece of code sooner or later
				for (Size j = 0; j < maps[i].size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					// DoubleReal rt_old = rt;
					transformations[i].apply(rt);
					maps[i][j].setRT(rt);
					// transformations[i].getPairs().push_back(TransformationDescription::PairVector::value_type(rt_old,rt));
				}
				// std::sort(transformations[i].getPairs().begin(),transformations[i].getPairs().end());
#endif
				setProgress(10*i+6);
			}
		}
		//set no transformation for reference map
		transformations[reference_map_index].setName("none");

		setProgress(10*maps.size());
		endProgress();
		return;
	}

	TransformationDescription MapAlignmentAlgorithmPoseClustering::calculateRegression_(Size const index_x_map, Size const index_y_map, ConsensusMap const& consensus_map, bool const symmetric_regression) const
	{
		// the result
		TransformationDescription lm;

		// coordinates of appropriate pairs will be stored here
		std::vector<double> vec_x;
		std::vector<double> vec_y;

		// search for pairs, optionally apply coordinate transformation, store coordinates in vec_x and vec_y
		FeatureHandle probe_x(index_x_map,0,Peak2D());
		FeatureHandle probe_y(index_y_map,0,Peak2D());
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
			if ( vec_x.size() == 1 ) // degenerate case, but we still can do something
			{
				slope = 1.;
				intercept = vec_y[0]-vec_x[0];
			}
			else // Oops, no pairs found at all!  Something went wrong outside.  Better throw an exception to indicate this fact.
			{
				throw Exception::Precondition
					(__FILE__,__LINE__,__PRETTY_FUNCTION__,
					 "not enough consensus features with feature handles from both maps"
					);
			}
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

} //namespace
