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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

//Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmApplyGivenTrafo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

// debugging yes/no
// #define V_MapAlignmentAlgorithm(a) std::cout << a << std::endl;
#define V_MapAlignmentAlgorithm(a)

using namespace std;

namespace OpenMS
{
	//register products here
	void MapAlignmentAlgorithm::registerChildren()
	{
		Factory<MapAlignmentAlgorithm>::registerProduct(
			MapAlignmentAlgorithmIdentification::getProductName(),
			&MapAlignmentAlgorithmIdentification::create);		

		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmPoseClustering::    getProductName(), &MapAlignmentAlgorithmPoseClustering::    create );
		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmSpectrumAlignment:: getProductName(), &MapAlignmentAlgorithmSpectrumAlignment:: create );
		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmApplyGivenTrafo::   getProductName(), &MapAlignmentAlgorithmApplyGivenTrafo::   create );
	}

	MapAlignmentAlgorithm::MapAlignmentAlgorithm()
		: DefaultParamHandler("MapAlignmentAlgorithm"),
			ProgressLogger()
	{
	}

	MapAlignmentAlgorithm::~MapAlignmentAlgorithm()
	{
	}

	void MapAlignmentAlgorithm::alignPeakMaps(vector< MSExperiment<> >&, vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}

	void MapAlignmentAlgorithm::alignFeatureMaps(vector< FeatureMap<> >&, vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);				
	}

	void MapAlignmentAlgorithm::alignConsensusMaps(vector<ConsensusMap>& cms, vector<TransformationDescription>& tf)
	{
    LOG_WARN << "MapAlignmentAlgorithm::alignConsensusMaps() does not support ConsensusMaps directly. Converting to FeatureMaps." << endl;

		vector< FeatureMap<> > maps_f;
    for (Size i=0; i<cms.size(); ++i)
    {
      FeatureMap<> fm;
      ConsensusMap::convert(cms[i], true, fm);
      maps_f.push_back(fm);
    }
    // call FeatureMap version of group()
    alignFeatureMaps(maps_f, tf);
    // apply transform
    transformConsensusMaps(cms, tf);
	}

	void MapAlignmentAlgorithm::alignPeptideIdentifications(vector< vector< PeptideIdentification > >&, vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);				
	}

	void MapAlignmentAlgorithm::setReference(Size reference_index,
																					 const String& reference_file)
	{
		if (reference_index || !reference_file.empty())
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This algorithm does not support a reference for the alignment.");
		}
	}

  void MapAlignmentAlgorithm::transformPeakMaps(
		vector< MSExperiment<> >& maps, 
		const vector<TransformationDescription>& given_trafos)
  {
    V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformPeakMaps()");

    if ( given_trafos.size() != maps.size() )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("MapAlignmentAlgorithm expects one given transformation (got: ") + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal");
    }

    for ( UInt index = 0; index < maps.size(); ++index )
    {
      MSExperiment<> & mse = maps[index];
      const TransformationDescription& td = given_trafos[index];
      transformSinglePeakMap(mse,td);
    }
    return;
  }

  void MapAlignmentAlgorithm::transformSinglePeakMap( MSExperiment<>& msexp, const TransformationDescription& trafo )
  {
    msexp.clearRanges();
    for ( MSExperiment<>::iterator mse_iter = msexp.begin(); mse_iter != msexp.end(); ++mse_iter )
    {
      DoubleReal rt = mse_iter->getRT();
      mse_iter->setRT(trafo.apply(rt));
    }
    msexp.updateRanges();
    return;
  }


  void MapAlignmentAlgorithm::transformFeatureMaps( vector< FeatureMap<> >& maps,
																										const vector<TransformationDescription>& given_trafos
		)
	{
		V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformFeatureMaps()");
		
		if ( given_trafos.size() != maps.size() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 String("MapAlignmentAlgorithm expects one given transformation (got: ")
				 + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
					);
		}

		for ( UInt index = 0; index < maps.size(); ++index )
		{
			FeatureMap<> & fm = maps[index];
			const TransformationDescription& td = given_trafos[index];
			transformSingleFeatureMap( fm, td );
		}

		return;
	}

	void MapAlignmentAlgorithm::transformSingleFeatureMap( FeatureMap<>& fmap, const TransformationDescription& trafo )
	{
		for ( vector<Feature>::iterator fmit = fmap.begin(); fmit != fmap.end(); ++fmit )
		{
			applyToFeature_(*fmit, trafo);
		}
		 
		// adapt RT values of unassigned peptides:
		if (!fmap.getUnassignedPeptideIdentifications().empty())
		{
			transformSinglePeptideIdentification(
				fmap.getUnassignedPeptideIdentifications(), trafo);
		}
		// if (trafo.getMaxRTErrorEstimate() > 0)
		// {
		//   LOG_WARN << "Maximal RT difference using alternative (linear) model was: " << trafo.getMaxRTErrorEstimate() << endl;
		// }
		return;
	}

	void MapAlignmentAlgorithm::applyToFeature_(
		Feature& feature, const TransformationDescription& trafo)
	{
		// V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::applyToFeature_()");

		applyToBaseFeature_(feature, trafo);

		// loop over all convex hulls
		vector<ConvexHull2D> & convex_hulls = feature.getConvexHulls();
		for ( vector<ConvexHull2D>::iterator chiter = convex_hulls.begin();
					chiter != convex_hulls.end(); ++chiter )
		{
			// transform all hull point positions within convex hull
			ConvexHull2D::PointArrayType points = chiter->getHullPoints();
			chiter->clear();
			for ( ConvexHull2D::PointArrayType::iterator points_iter = points.begin();
						points_iter != points.end();
						++points_iter
				)
			{
				DoubleReal rt = (*points_iter)[Feature::RT];
				(*points_iter)[Feature::RT] = trafo.apply(rt);
			}
			chiter->setHullPoints(points);
		}

		// recurse into subordinates
		for (vector<Feature>::iterator subiter = feature.getSubordinates().begin(); subiter != feature.getSubordinates().end(); ++subiter )
		{
			applyToFeature_(*subiter,trafo);
		}

		return;
	}


  void MapAlignmentAlgorithm::transformConsensusMaps(
		vector<ConsensusMap>& maps, 
		const vector<TransformationDescription>& given_trafos)
	{
		V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformConsensusMaps()");

		if ( given_trafos.size() != maps.size() )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("MapAlignmentAlgorithm expects one given transformation (got: ") + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal");
		}

		for (UInt index = 0; index < maps.size(); ++index)
		{
			ConsensusMap & cm = maps[index];
			const TransformationDescription& td = given_trafos[index];
			transformSingleConsensusMap(cm, td);
		}

		return;
	}

	void MapAlignmentAlgorithm::transformSingleConsensusMap(
		ConsensusMap& cmap, const TransformationDescription& trafo)
	{
		for (ConsensusMap::Iterator cmit = cmap.begin(); cmit != cmap.end(); 
				 ++cmit)
		{
			applyToConsensusFeature_(*cmit, trafo);
		}
    // if (trafo.getMaxRTErrorEstimate() > 0)
    // {
    //   LOG_WARN << "Maximal RT difference using alternative (linear) model was: " << trafo.getMaxRTErrorEstimate() << endl;
    // }

		// adapt RT values of unassigned peptides:
		if (!cmap.getUnassignedPeptideIdentifications().empty())
		{
			transformSinglePeptideIdentification(
				cmap.getUnassignedPeptideIdentifications(), trafo);
		}
		
		return;
	}


	void MapAlignmentAlgorithm::applyToBaseFeature_(
		BaseFeature& feature, const TransformationDescription& trafo)
	{
		// transform feature position:
		DoubleReal rt = feature.getRT();
		feature.setRT(trafo.apply(rt));

		// adapt RT values of annotated peptides:
		if (!feature.getPeptideIdentifications().empty())
		{
			transformSinglePeptideIdentification(feature.getPeptideIdentifications(),
																					 trafo);
		}
	}


	void MapAlignmentAlgorithm::applyToConsensusFeature_(
		ConsensusFeature& feature, const TransformationDescription& trafo )
	{
		// V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::applyToFeature_()");

		applyToBaseFeature_(feature, trafo);
		
		// apply to grouped features (feature handles):
		for (ConsensusFeature::HandleSetType::const_iterator it = 
					 feature.getFeatures().begin(); it != feature.getFeatures().end();
				 ++it)
		{
			DoubleReal rt = it->getRT();
			it->asMutable().setRT(trafo.apply(rt));
		}
		return;
	}


	void MapAlignmentAlgorithm::transformPeptideIdentifications(vector< vector< PeptideIdentification > >& maps, const vector<TransformationDescription>& given_trafos)
	{
		V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformPeptideIdentifications()");

		if ( given_trafos.size() != maps.size() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 String("MapAlignmentAlgorithm expects one given transformation (got: ")
				 + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
          );
		}

		for ( UInt map_index = 0; map_index < maps.size(); ++map_index )
		{
			V_MapAlignmentAlgorithm("map_index: " << map_index);
			const TransformationDescription& td = given_trafos[map_index];
			vector< PeptideIdentification >& pepids = maps[map_index];
			transformSinglePeptideIdentification(pepids,td);
		}
		return;
	}

	void MapAlignmentAlgorithm::transformSinglePeptideIdentification(vector< PeptideIdentification >& pepids, const TransformationDescription& trafo)
	{
		const UInt meta_index_RT = MetaInfo::registry().getIndex("RT");
		for ( UInt pepid_index = 0; pepid_index < pepids.size(); ++pepid_index )
		{
			V_MapAlignmentAlgorithm("pepid_index: " << pepid_index);
			PeptideIdentification & pepid = pepids[pepid_index];
			DataValue dv = pepid.getMetaValue(meta_index_RT);
			if (dv!=DataValue::EMPTY)
			{
				DoubleReal rt(dv);
				V_MapAlignmentAlgorithm("RT before: " << rt);
				rt = trafo.apply(rt);
				V_MapAlignmentAlgorithm("RT after: " << rt);
				pepid.setMetaValue(meta_index_RT, rt);
			}
		}
		return;
	}

	void MapAlignmentAlgorithm::fitModel(const String& model_type, const Param& params, vector<TransformationDescription>& trafos)
	{
		for (vector<TransformationDescription>::iterator it = trafos.begin();
				 it != trafos.end(); ++it)
		{
			it->fitModel(model_type, params);
		}
	}

	void MapAlignmentAlgorithm::getDefaultModel(String& model_type, Param& params)
	{
		model_type = "none";
		params.clear();
	}

} 
