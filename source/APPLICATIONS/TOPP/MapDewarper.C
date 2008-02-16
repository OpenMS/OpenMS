// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapDewarper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Grid.h>
#include <OpenMS/FORMAT/GridFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapDewarper MapDewarper
	
	@brief Dewarps a feature map by applying a transform to
	the coordinates of each feature.
	
	The dewarping is the last step in a map matching workflow. 
	The transform computed in the MapMatcher is now applied to the peak positions.
	Currently, we use a piecewise linear transform, but others can be implemented easily. This
	module simply applies this transform to the coordinates
	of each feature contained in the corresponding grid cells.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapDewarper
      : public TOPPBase
{
public:
  TOPPMapDewarper()
      : TOPPBase("MapDewarper","Dewarps a feature map by applying a transform to the coordinates")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
		registerInputFile_("feat","<file>","","the input FeatureXML file to be transformed");
		registerInputFile_("grid","<file>","","grid covering the map to be transformed");
    registerOutputFile_("out","<file>","","dewarped feature map");
  }

  ExitCodes main_(int , const char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

   	String gridfile = getStringOption_("grid");
    String features_file = getStringOption_("feat");
    String outfile = getStringOption_("out");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
		
    Grid the_grid;
    GridFile().load(gridfile,the_grid);

    FeatureXMLFile fmap_file;
    FeatureMap<> feature_map;
    fmap_file.load(features_file,feature_map);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    MapDewarper<> map_dewarper;
    map_dewarper.setMap(feature_map);
    map_dewarper.setGrid(the_grid);
    map_dewarper.dewarp();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    fmap_file.store(outfile,map_dewarper.getMap());

    return EXECUTION_OK;
  }
};


int main( int argc, const char** argv )
{
  TOPPMapDewarper tool;
  return tool.main(argc,argv);
}

/// @endcond
