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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LPWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LPWrapper* ptr = 0;
LPWrapper* null_ptr = 0;
START_SECTION(LPWrapper())
{
	ptr = new LPWrapper();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~LPWrapper())
{
	delete ptr;
}
END_SECTION

LPWrapper lp;
std::vector<DoubleReal> values(2,0.5);
std::vector<Int> indices;
indices.push_back(0);
indices.push_back(1);
START_SECTION((Size addColumn()))
{
  lp.addColumn();
  TEST_EQUAL(lp.getNumberOfColumns(),1);  
}
END_SECTION



START_SECTION((Size addRow(std::vector< Int > row_indices, std::vector< DoubleReal > row_values, String name)))
{
  lp.addColumn();
  lp.addRow(indices,values,String("row1"));
  TEST_EQUAL(lp.getNumberOfRows(),1);
  TEST_EQUAL(lp.getRowName(0),"row1");
}
END_SECTION


START_SECTION((Size addColumn(std::vector< Int > column_indices, std::vector< DoubleReal > column_values, String name)))
{
  lp.addRow(indices,values,String("row2"));
  lp.addColumn(indices,values,String("col3"));
  TEST_EQUAL(lp.getNumberOfColumns(),3);
  TEST_EQUAL(lp.getColumnName(2),"col3");
}
END_SECTION

START_SECTION((Size addRow(std::vector< Int > &row_indices, std::vector< DoubleReal > &row_values, String name, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.addRow(indices,values,String("row3"),0.2,1.2,LPWrapper::DOUBLE_BOUNDED_OR_FIXED);
  TEST_EQUAL(lp.getNumberOfRows(),3);
  TEST_EQUAL(lp.getRowName(2),"row3");
}
END_SECTION

START_SECTION((Size addColumn(std::vector< Int > &column_indices, std::vector< DoubleReal > &column_values, String name, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.addColumn(indices,values,String("col4"),0.2,1.2,LPWrapper::DOUBLE_BOUNDED_OR_FIXED);
  TEST_EQUAL(lp.getNumberOfColumns(),4);
  TEST_EQUAL(lp.getColumnName(3),"col4");
}
END_SECTION

START_SECTION((void setColumnName(Size index, String name)))
{
  lp.setColumnName(0,"col1");
  TEST_EQUAL(lp.getColumnName(0),"col1");  
}
END_SECTION

START_SECTION((String getColumnName(Size index)))
{
  TEST_EQUAL(lp.getColumnName(0),"col1"); 
}
END_SECTION

START_SECTION((String getRowName(Size index)))
{
  TEST_EQUAL(lp.getRowName(0),"row1"); 
}
END_SECTION

START_SECTION((Size getRowIndex(String name)))
{
  //  TEST_EQUAL(lp.getRowIndex("row1"),0); 
}
END_SECTION

START_SECTION((Size getColumnIndex(String name)))
{
  //TEST_EQUAL(lp.getColumnIndex("col1"),0); 
}
END_SECTION

START_SECTION((void setRowName(Size index, String name)))
{
  lp.setRowName(0,"new_row1");
  TEST_EQUAL(lp.getRowName(0),"new_row1"); 
}
END_SECTION

START_SECTION((void setColumnBounds(Size index, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.setColumnBounds(0,0.3,1.0,LPWrapper::DOUBLE_BOUNDED_OR_FIXED);
  TEST_EQUAL(lp.getColumnUpperBound(0),1.0);
  TEST_EQUAL(lp.getColumnLowerBound(0),0.3);
}
END_SECTION

START_SECTION((void setRowBounds(Size index, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.setRowBounds(0,-0.3,1.0,LPWrapper::DOUBLE_BOUNDED_OR_FIXED);
  TEST_EQUAL(lp.getRowUpperBound(0),1.0);
  TEST_EQUAL(lp.getRowLowerBound(0),-0.3);
}
END_SECTION

START_SECTION((void setColumnType(Size index, VariableType type)))
{
  lp.setColumnType(0,LPWrapper::INTEGER);
  TEST_EQUAL(lp.getColumnType(0),LPWrapper::INTEGER);
}
END_SECTION

START_SECTION((VariableType getColumnType(Size index)))
{
  lp.setColumnType(0,LPWrapper::BINARY);
  if (lp.getSolver()==LPWrapper::SOLVER_GLPK) TEST_EQUAL(lp.getColumnType(0),LPWrapper::BINARY)
  else TEST_EQUAL(lp.getColumnType(0),LPWrapper::INTEGER)
}
END_SECTION

START_SECTION((void setObjective(Size index, DoubleReal obj_value)))
{
  lp.setObjective(0,3.5);
  TEST_EQUAL(lp.getObjective(0),3.5);
}
END_SECTION

START_SECTION((DoubleReal getObjective(Size index)))
{
  lp.setObjective(1,2.5);
  TEST_EQUAL(lp.getObjective(1),2.5);
}
END_SECTION

START_SECTION((void setObjectiveSense(Sense sense)))
{
  lp.setObjectiveSense(LPWrapper::MIN);
  TEST_EQUAL(lp.getObjectiveSense(),LPWrapper::MIN)
}
END_SECTION

START_SECTION((Sense getObjectiveSense()))
{
  lp.setObjectiveSense(LPWrapper::MAX);
  TEST_EQUAL(lp.getObjectiveSense(),LPWrapper::MAX)
}
END_SECTION

START_SECTION((Size getNumberOfColumns()))
{
  TEST_EQUAL(lp.getNumberOfColumns(),4)
}
END_SECTION

START_SECTION((Size getNumberOfRows()))
{
  TEST_EQUAL(lp.getNumberOfRows(),3)
}
END_SECTION

START_SECTION((void readProblem(String filename, String format)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProblem(String filename, String format)))
{
  // TODO
}
END_SECTION

START_SECTION((Int solve(SolverParam &solver_param)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getStatus()))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getObjectiveValue()))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getColumnValue(Size index)))
{
  // TODO
}
END_SECTION

START_SECTION((void setSolver(const SOLVER s)))
{
  lp.setSolver(LPWrapper::SOLVER_GLPK);
  TEST_EQUAL(lp.getSolver(),LPWrapper::SOLVER_GLPK)
}
END_SECTION

START_SECTION((SOLVER getSolver() const ))
{
  TEST_EQUAL(lp.getSolver(),LPWrapper::SOLVER_GLPK)
}
END_SECTION

START_SECTION(([LPWrapper::SolverParam] SolverParam()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



