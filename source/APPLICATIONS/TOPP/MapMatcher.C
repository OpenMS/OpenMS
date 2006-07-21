// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>

#include <OpenMS/FORMAT/DGridFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>

#include <TOPPBase.h>

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
	
		
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES

class TOPPMapMatcher
            : public TOPPBase
{

public:
    TOPPMapMatcher()
            : TOPPBase("MapMatcher")
    {}

protected:

    void printToolUsage_()
    {
        cerr << endl
        << tool_name_ << " -- estimates a transformation for pairs of features in different LC/MS maps" << endl
        << endl
        << "Usage:" << endl
        << " " << tool_name_ << " [options]" << endl
        << endl
        << "Options are:" << endl
        << " -grid <file>   grid covering the map to be transformed (default read from INI file)" << endl
        << " -pairs <file>  feature pairs (default read from INI file)" << endl
        << " -q <float>  	 minimum quality of pairs considered (default read from INI file)" << endl
        << " -out <file>  	 output file (default read from INI file)" << endl
        << endl ;

    }

    void printToolHelpOpt_()
    {
        cerr << endl
        << tool_name_ << endl
        << endl
        << "INI options:" << endl
        << "  in <file>                  either feat or mzData (default read from INI file)" << endl
        << "  out <file>                output mzData file name" << endl
        << "  in_type <file_type>  either feat or mzData (default read from INI file)" << endl
        << endl
        << "INI File example section:" << endl
        << "<ITEM name=\"pairs\" value=\"04111717_pairs.xml\" type=\"string\"/>" << endl
        << "<ITEM name=\"grid\" value=\"the_grid.xml\" type=\"string\"/>" << endl
        << "<ITEM name=\"min_quality\" value=\"0.5\" type=\"float\"/>" << endl
        << "<ITEM name=\"out\" value=\"grid_wtransform.xml\" type=\"string\"/>" << endl;
    }

    void setOptionsAndFlags_()
    {
        options_["-grid"] = "grid";
        options_["-pairs"] = "pairs";
        options_["-out"] = "out";
        options_["-q"] = "min_quality";
    }

    ExitCodes main_(int , char**)
    {
        //-----------------------------------------------------------
        // parsing parameters
        //-------------------------------------------------------------
        //File names
        String gridfile, pairsfile, outfile;

        gridfile = getParamAsString_("grid");
        writeDebug_(String("Grid file: ") + gridfile, 1);

        pairsfile = getParamAsString_("pairs");
        writeDebug_(String("Pairs file: ") + pairsfile, 1);

        outfile = getParamAsString_("out");
        writeDebug_(String("Output file: ") + outfile, 1);

        //parameters
        double min_quality = 0;

        //length of the structuring element
        String qualstr = getParamAsString_("min_quality");
        writeDebug_(String("min_quality") + qualstr, 1);

        try
        {
            //resampling
            if (qualstr != "")
            {
                min_quality = qualstr.toFloat();
            }
        }
        catch(Exception::ConversionError& e)
        {
            writeLog_(String("Invalid value for the minimum quality '") + qualstr  + "' given. Aborting!");
            printUsage_();
            return ILLEGAL_PARAMETERS;
        }

        //-------------------------------------------------------------
        // reading input
        //-------------------------------------------------------------

        DGridFile grid_file;
        DGrid<2> the_grid;
        grid_file.load(gridfile,the_grid);

        DFeaturePairsFile pairs_file;
        DFeaturePairVector<2> pairs_vector;
        pairs_file.load(pairsfile,pairs_vector);

        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------
        DMapMatcherRegression<2> map_matcher;
        map_matcher.setFeaturePairs(pairs_vector);
        map_matcher.setGrid(the_grid);
        map_matcher.setMinQuality(min_quality);

        map_matcher.estimateTransform();

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        DGrid<2> grid_with_transform;
        grid_with_transform = map_matcher.getGrid();

        grid_file.store(outfile,grid_with_transform);

        return OK;

    }
};

/// @endcond

int main( int argc, char ** argv )
{
    TOPPMapMatcher tool;
    return tool.main(argc,argv);
}





