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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DMapDewarper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
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
	@page MapDewarper MapDewarper
	
	@brief Dewarps a feature map by applying a transform to
	the coordinates of each feature.
	
	The dewarping is the last and optional step in a map
	matching workflow. The transform was computed in the map matching
	step of the workflow. Currently, we use a piecewise
	linear transform, but others can be implemented easily. This
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
		registerStringOption_("feat","<file>","","the feature map to be transformed");
		registerStringOption_("grid","<file>","","grid covering the map to be transformed");
    registerStringOption_("out","<file>","","dewarped feature map");
  }

  ExitCodes main_(int , char**)
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
		
    DGrid<2> the_grid;
    DGridFile().load(gridfile,the_grid);

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

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    fmap_file.store(outfile,map_dewarper.getMap());

    return EXECUTION_OK;
  }
};


int main( int argc, char ** argv )
{
  TOPPMapDewarper tool;
  return tool.main(argc,argv);
}

/// @endcond
