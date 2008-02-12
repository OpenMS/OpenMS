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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include<OpenMS/FORMAT/GridFile.h>

#include<OpenMS/ANALYSIS/MAPMATCHING/Grid.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/GridCell.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

#include<vector>

///////////////////////////

START_TEST(GridFile_test, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

GridFile* ptr = 0;
CHECK((GridFile()))
	ptr = new GridFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~GridFile()))
	delete ptr;
RESULT

CHECK((void load(String filename, Grid &grid) throw (Exception::FileNotFound, Exception::ParseError)))
	Grid grid;
  GridFile().load("data/GridFile.xml",grid);
	
	//test size
  TEST_EQUAL(grid.size(),1);
  //test bounds
	TEST_REAL_EQUAL(grid.back().minX(),1);
	TEST_REAL_EQUAL(grid.back().minY(),2);
	TEST_REAL_EQUAL(grid.back().maxX(),10);
	TEST_REAL_EQUAL(grid.back().maxY(),11);
	//test mappings
	TEST_EQUAL(grid.back().getMappings().size(),2);
	LinearMapping* mapping = dynamic_cast<LinearMapping*>(grid.back().getMappings()[0]);
	TEST_NOT_EQUAL(mapping, 0);
	mapping = dynamic_cast<LinearMapping*>(grid.back().getMappings()[1]);
	TEST_NOT_EQUAL(mapping, 0);
RESULT

CHECK((void store(String filename, const Grid &grid) const  throw (Exception::UnableToCreateFile)))
	
	std::string tmp_filename;
  Grid grid;
	GridFile gfile;
  
  NEW_TMP_FILE(tmp_filename);
  gfile.load("data/GridFile.xml",grid);
	gfile.store(tmp_filename,grid);
	
	TEST_FILE(tmp_filename.c_str(),"data/GridFile.xml");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
