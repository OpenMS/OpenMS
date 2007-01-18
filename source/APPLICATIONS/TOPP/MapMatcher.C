// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// Initial version by Ole Schulz-Trieglaff
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/FORMAT/DGridFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapMatcher MapMatcher

	@brief Computes a transformation for a list of feature pairs.

	This is the second step in the map matching workflow. This application
	takes a list of feature pairs as computed by the FeatureMatcher and
	a grid (partially) covering the LC/MS map. For each grid cell, a
	transformation is computed that maps the feature partners on each
	other. Currently, this transformation is linear.

	The output of this application is the list of grid cells with the
	estimated transformation.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapMatcher
	: public TOPPBase
{

 public:
	TOPPMapMatcher()
		: TOPPBase("MapMatcher", "estimate a transformation to map a list of pairs of features in different LC/MS maps onto each other")
	{}

 protected:

	void registerOptionsAndFlags_()
	{
		registerStringOption_("pairs","<file>","","input feature pairs file");
		registerStringOption_("grid","<file>","","input grid file");
		registerStringOption_("out","<file>","","output grid file");
		registerDoubleOption_("min_quality","<double>",0,"minimum quality of pairs considered",false);
	}

	ExitCodes main_(int , char**)
	{
		//-----------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------

		String grid_filename = getStringOption_("grid");
		inputFileReadable_(grid_filename);

		String pairs_filename = getStringOption_("pairs");
		inputFileReadable_(pairs_filename);

		String out_filename = getStringOption_("out");
		outputFileWritable_(out_filename);

		double min_quality = getDoubleOption_("min_quality");

		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------

		DGridFile grid_file;
		DGrid<2> grid;
		writeLog_("Reading grid file " + grid_filename );
		grid_file.load(grid_filename, grid);

		DFeaturePairsFile pairs_file;
		DFeaturePairVector<2> pairs_vector;
		writeLog_("Reading pairs file " + pairs_filename );
		pairs_file.load(pairs_filename, pairs_vector);

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		DMapMatcherRegression<> map_matcher;
		map_matcher.setFeaturePairs(pairs_vector);
		map_matcher.setGrid(grid);
		map_matcher.setMinQuality(min_quality);

		// action!
		map_matcher.estimateTransform();

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------

		DGrid<2> const & grid_with_transform = map_matcher.getGrid();

		grid_file.store(out_filename, grid_with_transform);

		return EXECUTION_OK;

	}
};

/// @endcond

int main( int argc, char ** argv )
{
	TOPPMapMatcher tool;
	return tool.main(argc,argv);
}
