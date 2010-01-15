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
// $Authors: $
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
