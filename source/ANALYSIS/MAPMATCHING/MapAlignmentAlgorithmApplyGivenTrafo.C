// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmApplyGivenTrafo.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

// debugging yes/no
// #define V_MapAlignmentAlgorithmApplyGivenTrafo(a) std::cout << a << std::endl;
#define V_MapAlignmentAlgorithmApplyGivenTrafo(a)

namespace OpenMS
{

	MapAlignmentAlgorithmApplyGivenTrafo::MapAlignmentAlgorithmApplyGivenTrafo()
		: MapAlignmentAlgorithm()
	{
		setName("MapAlignmentAlgorithmApplyGivenTrafo");
	
		defaults_.setValue("transformations",StringList(),"file names of transformations to be used (TrafoXML format)");
		defaults_.setValue("transformations_path","","optional path prepended to all transformations");
		// defaults_.setValue("inverse","","do the inverse transformation"); // TODO would not this be nice.
		defaultsToParam_();
		
		return;
	}

	MapAlignmentAlgorithmApplyGivenTrafo::~MapAlignmentAlgorithmApplyGivenTrafo()
	{
		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::setGivenTrafos(const std::vector<TransformationDescription>& given_trafos)
	{
		given_trafos_ = given_trafos;
		return;
	}

	std::vector<TransformationDescription>& MapAlignmentAlgorithmApplyGivenTrafo::getGivenTrafos()
	{
		return given_trafos_;
	}
	
	const std::vector<TransformationDescription>& MapAlignmentAlgorithmApplyGivenTrafo::getGivenTrafos() const
	{
		return given_trafos_;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::alignPeakMaps( std::vector< MSExperiment<> >& maps,
																														std::vector<TransformationDescription>& transformations
																													)
	{
		std::cout << "Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::alignPeakMaps()" << std::endl;

		if ( !transformations.empty() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 "MapAlignmentAlgorithmApplyGivenTrafo does not output transformations, the list must be empty"
				);
		}
		
		readGivenTrafos();

		transformPeakMaps( maps, given_trafos_ );

		return;
	}
		
	void MapAlignmentAlgorithmApplyGivenTrafo::transformPeakMaps( std::vector< MSExperiment<> >& maps,
																																const  std::vector<TransformationDescription>& given_trafos
																															)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::transformPeakMaps()");

		if ( given_trafos.size() != maps.size() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 String("MapAlignmentAlgorithmApplyGivenTrafo expects one given transformation (got: ")
				 + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
				);
		}

		for ( UInt index = 0; index < maps.size(); ++index )
		{
			MSExperiment<> & mse = maps[index];
			mse.clearRanges();
			const TransformationDescription& td = given_trafos[index];
			td.init_();
			for ( MSExperiment<>::iterator mse_iter = mse.begin(); mse_iter != mse.end(); ++mse_iter )
			{
				DoubleReal rt = mse_iter->getRT();
				(*td.trafo_)(rt);
				mse_iter->setRT(rt);
			}
			mse.updateRanges();
		}
		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::alignFeatureMaps( std::vector< FeatureMap<> >& maps,
																															 std::vector<TransformationDescription>& transformations
																														 )
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::alignFeatureMaps()");

		if ( !transformations.empty() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 "MapAlignmentAlgorithmApplyGivenTrafo does not output transformations, the list must be empty"
				);
		}
		
		readGivenTrafos();

		transformFeatureMaps( maps, given_trafos_ );

		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::transformFeatureMaps( std::vector< FeatureMap<> >& maps,
																																	 const std::vector<TransformationDescription>& given_trafos
																																 )
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::transformFeatureMaps()");

		if ( given_trafos.size() != maps.size() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 String("MapAlignmentAlgorithmApplyGivenTrafo expects one given transformation (got: ")
				 + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
				);
		}

		for ( UInt index = 0; index < maps.size(); ++index )
		{
			FeatureMap<> & fm = maps[index];
			const TransformationDescription& td = given_trafos[index];
			td.init_();
			for ( std::vector<Feature>::iterator fmit = fm.begin(); fmit != fm.end(); ++fmit )
			{
				applyToFeature_(fmit,*td.trafo_);
			}
		}

		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::applyToFeature_( const std::vector<Feature>::iterator &iter,
																															TransformationDescription::Trafo_ const& trafo
																														)
	{
		// V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::applyToFeature_()");

		// transform feature position
		DoubleReal rt = iter->getRT();
		trafo(rt);
		iter->setRT(rt);

		// loop over all convex hulls
		std::vector<ConvexHull2D> & convex_hulls = iter->getConvexHulls();
		for ( std::vector<ConvexHull2D>::iterator chiter = convex_hulls.begin();
					chiter!= convex_hulls.end();
					++chiter
				)
		{
			// transform all hull point positions within convex hull
			ConvexHull2D::PointArrayType & points = const_cast<ConvexHull2D::PointArrayType&>(chiter->getPoints());
			for ( ConvexHull2D::PointArrayType::iterator points_iter = points.begin();
						points_iter != points.end();
						++points_iter
					)
			{
				DoubleReal rt = (*points_iter)[Feature::RT];
				trafo(rt);
				(*points_iter)[Feature::RT] = rt;
			}
		}

		// recurse into subordinates
		for ( std::vector<Feature>::iterator subiter = iter->getSubordinates().begin();
					subiter != iter->getSubordinates().end();
					++subiter )
		{
			applyToFeature_(subiter,trafo);
		}

		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::alignPeptideIdentifications( std::vector< std::vector< PeptideIdentification > >& maps,
																																					std::vector<TransformationDescription>& transformations
																																				)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::alignPeptideIdentifications()");
		
		if ( !transformations.empty() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 "MapAlignmentAlgorithmApplyGivenTrafo does not output transformations, the list must be empty"
				);
		}
		
		readGivenTrafos();

		transformPeptideIdentifications( maps, given_trafos_ );

		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::transformPeptideIdentifications( std::vector< std::vector< PeptideIdentification > >& maps,
																																							const std::vector<TransformationDescription>& given_trafos
																																						)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::transformPeptideIdentifications()");

		if ( given_trafos.size() != maps.size() )
		{
			throw Exception::IllegalArgument
				(__FILE__, __LINE__, __PRETTY_FUNCTION__,
				 String("MapAlignmentAlgorithmApplyGivenTrafo expects one given transformation (got: ")
				 + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
				);
		}

		const UInt meta_index_RT = MetaInfo::registry().getIndex("RT");

		for ( UInt map_index = 0; map_index < maps.size(); ++map_index )
		{
			V_MapAlignmentAlgorithmApplyGivenTrafo("map_index: " << map_index);
			const TransformationDescription& td = given_trafos[map_index];
			td.init_();
			for ( UInt pepid_index = 0; pepid_index < maps[map_index].size(); ++pepid_index )
			{
				V_MapAlignmentAlgorithmApplyGivenTrafo("pepid_index: " << pepid_index);
				PeptideIdentification & pepid = maps[map_index][pepid_index];
				DataValue dv = pepid.getMetaValue(meta_index_RT);
				if (dv!=DataValue::EMPTY)
				{
					DoubleReal rt(dv);
					V_MapAlignmentAlgorithmApplyGivenTrafo("RT before: " << rt);
					(*td.trafo_)(rt);
					pepid.setMetaValue(meta_index_RT,rt);
					V_MapAlignmentAlgorithmApplyGivenTrafo("RT after: " << rt);
				}
			}
		}

		return;
	}

	void MapAlignmentAlgorithmApplyGivenTrafo::readGivenTrafos()
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::readGivenTrafos()");

		StringList given_trafo_files = param_.getValue("transformations");
		String transformations_path = param_.getValue("transformations_path");
		TransformationXMLFile trafo_xml;
		for ( StringList::const_iterator slcit = given_trafo_files.begin();
					slcit != given_trafo_files.end();
					++slcit
				)
		{
			given_trafos_.push_back(TransformationDescription());
			trafo_xml.load( transformations_path + *slcit, given_trafos_.back() );
			// std::cout << "I just loaded this trafo:\n" << given_trafos.back() << "\n" << std::endl;
		}
		return;
	}


} //namespace 
