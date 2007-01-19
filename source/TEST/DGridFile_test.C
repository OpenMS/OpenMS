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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include<OpenMS/FORMAT/DGridFile.h>

#include<OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>

#include<vector>

///////////////////////////

START_TEST(DGridFile_test, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DGridFile* ptr = 0;
CHECK((DGridFile()))
	ptr = new DGridFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~DGridFile()))
	delete ptr;
RESULT

CHECK((template<Size D> void load(String filename, DGrid<D>& grid) throw(Exception::FileNotFound, Exception::ParseError)))
	PRECISION(0.01)
	
	DGrid<2> grid;
	DGridFile gfile;
		   
  gfile.load("data/DGridFile.xml",grid);
  DGridCell<2> cell = grid.back();
  	
	TEST_EQUAL(cell.minX(),0);
	TEST_EQUAL(cell.minY(),0);
	TEST_EQUAL(cell.maxX(),10);
	TEST_EQUAL(cell.maxY(),10);
	
	DGridCell<2>::MappingVector mappings = cell.getMappings();
	
	TEST_EQUAL(mappings.size(),2);
	
//	DLinearMapping<1>* mapp1 = static_cast<DLinearMapping<1>*>(mappings.back());
//	DLinearMapping<1>* mapp2 = static_cast<DLinearMapping<1>*>(mappings.back());
	
	
RESULT

CHECK((template<Size D> void store(String filename, const DGrid<D>& grid) const throw(Exception::UnableToCreateFile)))
	
	std::string tmp_filename;
  DGrid<2> grid;
	DGridFile gfile;
  
  NEW_TMP_FILE(tmp_filename);
  gfile.load("data/DGridFile.xml",grid);
	gfile.store(tmp_filename,grid);
	
	TEST_FILE(tmp_filename.c_str(),"data/DGridFile.xml");
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
