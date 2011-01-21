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
	
		defaults_.setValue("transformations", StringList(), "file names of transformations to be used (TrafoXML format)");
		defaults_.setValue("transformations_path", "", "optional path prepended to all transformations");
		defaults_.setValue("invert", "false", "compute and apply the inverse transformation");
		defaults_.setValidStrings("invert", StringList::create("true,false"));
		defaultsToParam_();
	}

	MapAlignmentAlgorithmApplyGivenTrafo::~MapAlignmentAlgorithmApplyGivenTrafo()
	{
	}


	void MapAlignmentAlgorithmApplyGivenTrafo::alignPeakMaps(std::vector< MSExperiment<> >& /* maps */, std::vector<TransformationDescription>& transformations)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::alignPeakMaps()");
		
		readGivenTrafos(transformations);
	}


	void MapAlignmentAlgorithmApplyGivenTrafo::alignFeatureMaps(std::vector< FeatureMap<> >& /* maps */, std::vector<TransformationDescription>& transformations)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::alignFeatureMaps()");

		readGivenTrafos(transformations);
	}


	void MapAlignmentAlgorithmApplyGivenTrafo::alignPeptideIdentifications(std::vector< std::vector< PeptideIdentification > >& /* maps */, std::vector<TransformationDescription>& transformations)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there.  This is MapAlignmentAlgorithmApplyGivenTrafo::alignPeptideIdentifications()");
		
		readGivenTrafos(transformations);
	}


	void MapAlignmentAlgorithmApplyGivenTrafo::readGivenTrafos(std::vector<TransformationDescription>& transformations)
	{
		V_MapAlignmentAlgorithmApplyGivenTrafo("Hi out there. This is MapAlignmentAlgorithmApplyGivenTrafo::readGivenTrafos()");

		StringList given_trafo_files = param_.getValue("transformations");
		String transformations_path = param_.getValue("transformations_path");
		TransformationXMLFile trafo_xml;
		bool invert = (String(param_.getValue("invert")) == "true");

		transformations.resize(given_trafo_files.size());
		for (Size i = 0; i < given_trafo_files.size(); ++i)
		{		
			trafo_xml.load(transformations_path + given_trafo_files[i], 
										 transformations[i]);
			if (invert)
			{
				transformations[i].invert();
			}
		}
	}


	void MapAlignmentAlgorithmApplyGivenTrafo::getDefaultModel(String& model_type, Param& params)
	{
		model_type = "none";
		params.clear();
	}

} //namespace 
