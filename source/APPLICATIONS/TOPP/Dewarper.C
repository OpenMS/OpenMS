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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DMapDewarper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>

#include <OpenMS/FORMAT/DGridFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/KERNEL/DFeatureMap.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page Dewarper Dewarper
	
	@brief Dewarps a feature map by applying a transform to
	the coordinates of each feature.
	
	The dewarping is the last and optional step in a map
	matching workflow. The transform was computed in the map matching
	step of the workflow. Currently, we use a piecewise
	linear transform, but others can be implemented easily. This
	module simply applies this transform to the coordinates
	of each feature contained in the corresponding grid cells.
		
	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDewarper
      : public TOPPBase
{
public:
  TOPPDewarper()
      : TOPPBase("Dewarper")
  {}

protected:
  void printToolUsage_() const
  {
    cerr << endl
    << getToolName() << " -- dewarps a feature map" << endl
    << "Version: " << VersionInfo::getVersion() << endl
    << endl
    << "Usage:" << endl
    << " " << getToolName() << " [options]" << endl
    << endl
    << "Options are:" << endl
    << "  -grid <file>   grid covering the map to be transformed" << endl
    << "  -feat <file>   feature pairs" << endl
    << "  -out <file>    dewarped feature map" << endl
    << endl;
  }

  void printToolHelpOpt_() const
  {
    cerr << endl
    << getToolName() << endl
    << endl
    << "INI options:" << endl
    << "grid  grid covering the map to be transformed" << endl
    << "  feat  feature pairs" << endl
    << "  out   dewarped feature map" << endl
    << endl
    << "INI File example section:" << endl
    << "  <ITEM name=\"grid\" value=\"grid.xml\" type=\"string\"/>" << endl
    << "  <ITEM name=\"feat\" value=\"input.feat\" type=\"string\"/>" << endl
    << "  <ITEM name=\"out\" value=\"output.feat\" type=\"string\"/>" << endl;
  }

  void setOptionsAndFlags_()
  {
    options_["-out"] = "out";
    options_["-grid"] = "grid";
    options_["-feat"] = "feat";
  }

  ExitCodes main_(int , char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String gridfile, features_file, outfile;
	   
	gridfile = getParamAsString_("grid");
    writeDebug_(String("Grid file: ") + gridfile,1);

    features_file = getParamAsString_("feat");
    writeDebug_(String("Feature file: ")+ features_file,1);

    //determine output file name
    outfile = getParamAsString_("out");
    writeDebug_(String("Output file: ")+ outfile,1);


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    DGridFile grid_file;
    DGrid<2> the_grid;
    grid_file.load(gridfile,the_grid);

    DFeatureMapFile fmap_file;
    DFeatureMap<2> feature_map;
    fmap_file.load(features_file,feature_map);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    DMapDewarper<> map_dewarper;
    map_dewarper.setMap(feature_map);
    map_dewarper.setGrid(the_grid);

    map_dewarper.dewarp();

    DFeatureMap<2> dewarped_features;
    dewarped_features = map_dewarper.getMap();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    fmap_file.store(outfile,dewarped_features);

    return EXECUTION_OK;
  }
};


int main( int argc, char ** argv )
{
  TOPPDewarper tool;
  return tool.main(argc,argv);
}

/// @endcond
