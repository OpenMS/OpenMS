// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#if COINOR_SOLVER==1
  #include "coin/CoinModel.hpp"
#endif
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
//lp.setSolver(LPWrapper::SOLVER_GLPK);
std::vector<DoubleReal> values(2,0.5);
std::vector<Int> indices;
indices.push_back(0);
indices.push_back(1);
START_SECTION((Int addColumn()))
{
  lp.addColumn();
  TEST_EQUAL(lp.getNumberOfColumns(),1);  
}
END_SECTION



START_SECTION((Int addRow(std::vector< Int > row_indices, std::vector< DoubleReal > row_values, const String &name)))
{
  lp.addColumn();
  lp.addRow(indices,values,String("row1"));
  TEST_EQUAL(lp.getNumberOfRows(),1);
  TEST_EQUAL(lp.getRowName(0),"row1");
}
END_SECTION


START_SECTION((Int addColumn(std::vector< Int > column_indices, std::vector< DoubleReal > column_values, const String &name)))
{
  lp.addRow(indices,values,String("row2"));
  lp.addColumn(indices,values,String("col3"));
  TEST_EQUAL(lp.getNumberOfColumns(),3);
  TEST_EQUAL(lp.getColumnName(2),"col3");
}
END_SECTION

START_SECTION((Int addRow(std::vector< Int > &row_indices, std::vector< DoubleReal > &row_values, const String &name, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.addRow(indices,values,String("row3"),0.2,1.2,LPWrapper::DOUBLE_BOUNDED);
  TEST_EQUAL(lp.getNumberOfRows(),3);
  TEST_EQUAL(lp.getRowName(2),"row3");
}
END_SECTION

START_SECTION((Int addColumn(std::vector< Int > &column_indices, std::vector< DoubleReal > &column_values, const String &name, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.addColumn(indices,values,String("col4"),0.2,1.2,LPWrapper::DOUBLE_BOUNDED);
  TEST_EQUAL(lp.getNumberOfColumns(),4);
  TEST_EQUAL(lp.getColumnName(3),"col4");
}
END_SECTION

START_SECTION((void setColumnName(Int index, const String &name)))
{
  lp.setColumnName(0,"col1");
  TEST_EQUAL(lp.getColumnName(0),"col1");  
}
END_SECTION

START_SECTION((String getColumnName(Int index)))
{
  TEST_EQUAL(lp.getColumnName(0),"col1"); 
}
END_SECTION

START_SECTION((String getRowName(Int index)))
{
  TEST_EQUAL(lp.getRowName(0),"row1"); 
}
END_SECTION

START_SECTION((Int getRowIndex(const String &name)))
{
  TEST_EQUAL(lp.getRowIndex("row1"),0); 
}
END_SECTION

START_SECTION((Int getColumnIndex(const String &name)))
{
  TEST_EQUAL(lp.getColumnIndex("col1"),0); 
}
END_SECTION

START_SECTION((void setRowName(Int index, const String &name)))
{
  lp.setRowName(0,"new_row1");
  TEST_EQUAL(lp.getRowName(0),"new_row1"); 
}
END_SECTION

START_SECTION((void setColumnBounds(Int index, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.setColumnBounds(0,0.3,1.0,LPWrapper::DOUBLE_BOUNDED);
  TEST_EQUAL(lp.getColumnUpperBound(0),1.0);
  TEST_EQUAL(lp.getColumnLowerBound(0),0.3);
}
END_SECTION

START_SECTION((void setRowBounds(Int index, DoubleReal lower_bound, DoubleReal upper_bound, Type type)))
{
  lp.setRowBounds(0,-0.3,1.0,LPWrapper::DOUBLE_BOUNDED);
  TEST_EQUAL(lp.getRowUpperBound(0),1.0);
  TEST_EQUAL(lp.getRowLowerBound(0),-0.3);
}
END_SECTION

START_SECTION((void setColumnType(Int index, VariableType type)))
{
  lp.setColumnType(0,LPWrapper::INTEGER);
  TEST_EQUAL(lp.getColumnType(0),LPWrapper::INTEGER);
}
END_SECTION

START_SECTION((VariableType getColumnType(Int index)))
{
  lp.setColumnType(1,LPWrapper::BINARY);
  if (lp.getSolver()==LPWrapper::SOLVER_GLPK) TEST_EQUAL(lp.getColumnType(1),LPWrapper::BINARY)
  else TEST_EQUAL(lp.getColumnType(1),LPWrapper::INTEGER)
}
END_SECTION

START_SECTION((void setObjective(Int index, DoubleReal obj_value)))
{
  lp.setObjective(0,3.5);
  TEST_EQUAL(lp.getObjective(0),3.5);
}
END_SECTION

START_SECTION((DoubleReal getObjective(Int index)))
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

START_SECTION((Int getNumberOfColumns()))
{
  TEST_EQUAL(lp.getNumberOfColumns(),4)
}
END_SECTION

START_SECTION((Int getNumberOfRows()))
{
  TEST_EQUAL(lp.getNumberOfRows(),3)
}
END_SECTION

START_SECTION((DoubleReal getColumnUpperBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getColumnUpperBound(0),1.0)
}
END_SECTION

START_SECTION((void deleteRow(Int index)))
{
  lp.deleteRow(2);
  if(lp.getSolver() == LPWrapper::SOLVER_GLPK)
    {
      TEST_EQUAL(lp.getNumberOfRows(),2)
    }
#if COINOR_SOLVER==1
  else
    {
      // CoinOr doesn't delete the column, but sets all entries to zero and deletes the bounds, names, objective coeff etc.
      TEST_REAL_SIMILAR(lp.getObjective(2),0.)
      TEST_REAL_SIMILAR(lp.getColumnLowerBound(2),-COIN_DBL_MAX)
      TEST_REAL_SIMILAR(lp.getColumnUpperBound(2),COIN_DBL_MAX)  
    }
#endif
}
END_SECTION

START_SECTION((DoubleReal getColumnLowerBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getColumnLowerBound(0),0.3)
}
END_SECTION

START_SECTION((DoubleReal getRowUpperBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getRowUpperBound(0),1.0)
}
END_SECTION

START_SECTION((DoubleReal getRowLowerBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getRowLowerBound(0),-0.3)  
}
END_SECTION


START_SECTION((void setElement(Int row_index, Int column_index, DoubleReal value)))
{
  lp.setElement(1,2,0.5);
  TEST_REAL_SIMILAR(lp.getElement(1,2),0.5)  
}
END_SECTION


START_SECTION((DoubleReal getElement(Int row_index, Int column_index)))
{
  lp.setElement(0,2,0.1);
  TEST_REAL_SIMILAR(lp.getElement(0,2),0.1)
}
END_SECTION


START_SECTION((void readProblem(String filename, String format)))
{
  if(lp.getSolver() == LPWrapper::SOLVER_GLPK)
    {
      lp.readProblem(OPENMS_GET_TEST_DATA_PATH("LPWrapper_test.lp"),"LP");
      TEST_EQUAL(lp.getNumberOfColumns(),2)
      TEST_EQUAL(lp.getNumberOfRows(),3)
      TEST_EQUAL(lp.getColumnType(0),LPWrapper::INTEGER)
      TEST_EQUAL(lp.getColumnType(1),LPWrapper::INTEGER)
      TEST_EQUAL(lp.getObjective(0),1)
      TEST_EQUAL(lp.getObjective(1),0)
      TEST_EQUAL(lp.getRowUpperBound(0),0)
      TEST_EQUAL(lp.getRowUpperBound(1),12)
      TEST_EQUAL(lp.getRowUpperBound(2),12)
      TEST_EQUAL(lp.getElement(0,0),1)
      TEST_EQUAL(lp.getElement(0,1),-1)
      TEST_EQUAL(lp.getElement(1,0),2)
      TEST_EQUAL(lp.getElement(1,1),3)
      TEST_EQUAL(lp.getElement(2,0),3)
      TEST_EQUAL(lp.getElement(2,1),2) 
    }
#if COINOR_SOLVER==1
  else  if (lp.getSolver()==LPWrapper::SOLVER_COINOR)
    {
      TEST_EXCEPTION(Exception::NotImplemented, lp.readProblem("/bla/bluff/blblb/sdfhsdjf/test.txt","LP"))
    }
#endif
}
END_SECTION

START_SECTION((void writeProblem(const String &filename, const WriteFormat format) const ))
{
  if(lp.getSolver() == LPWrapper::SOLVER_GLPK)
    {
      String tmp_filename;
      NEW_TMP_FILE(tmp_filename);
      lp.writeProblem(tmp_filename,LPWrapper::FORMAT_LP);
      LPWrapper lp2;
      lp2.setSolver(LPWrapper::SOLVER_GLPK);
      lp2.readProblem(tmp_filename,"LP");
      TEST_EQUAL(lp.getNumberOfColumns(),2)
      TEST_EQUAL(lp.getNumberOfRows(),3)
      TEST_EQUAL(lp.getColumnType(0),LPWrapper::INTEGER)
      TEST_EQUAL(lp.getColumnType(1),LPWrapper::INTEGER)
      TEST_EQUAL(lp.getObjective(0),1)
      TEST_EQUAL(lp.getObjective(1),0)
      TEST_EQUAL(lp.getRowUpperBound(0),0)
      TEST_EQUAL(lp.getRowUpperBound(1),12)
      TEST_EQUAL(lp.getRowUpperBound(2),12)
      TEST_EQUAL(lp.getElement(0,0),1)
      TEST_EQUAL(lp.getElement(0,1),-1)
      TEST_EQUAL(lp.getElement(1,0),2)
      TEST_EQUAL(lp.getElement(1,1),3)
      TEST_EQUAL(lp.getElement(2,0),3)
      TEST_EQUAL(lp.getElement(2,1),2)
    }
#if COINOR_SOLVER==1
  else  if (lp.getSolver()==LPWrapper::SOLVER_COINOR)
  {
    TEST_EXCEPTION(Exception::NotImplemented, lp.writeProblem("/bla/bluff/blblb/sdfhsdjf/test.txt",LPWrapper::FORMAT_LP))
  }
#endif
}
END_SECTION

START_SECTION((Int solve(SolverParam &solver_param, const Size verbose_level=0)))
{
#if COINOR_SOLVER==1
  if(lp.getSolver() ==LPWrapper::SOLVER_COINOR)   lp.readProblem(OPENMS_GET_TEST_DATA_PATH("LPWrapper_test.mps"),"MPS");
  lp.setObjectiveSense(LPWrapper::MAX);
#endif
  LPWrapper::SolverParam param;
  lp.solve(param);
}
END_SECTION

START_SECTION((SolverStatus getStatus()))
{
  if(lp.getSolver() == LPWrapper::SOLVER_GLPK)
    {
      TEST_EQUAL(lp.getStatus(),LPWrapper::OPTIMAL)
    }
#if COINOR_SOLVER==1
  else
  {
    TEST_EQUAL(lp.getStatus(),LPWrapper::UNDEFINED)
  }
#endif

}
END_SECTION

START_SECTION((DoubleReal getObjectiveValue()))
{
  TEST_REAL_SIMILAR(lp.getObjectiveValue(),2)
}
END_SECTION

START_SECTION((DoubleReal getColumnValue(Int index)))
{
  TEST_REAL_SIMILAR(lp.getColumnValue(0),2)
  TEST_REAL_SIMILAR(lp.getColumnValue(1),2)
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
  LPWrapper::SolverParam* sptr = new LPWrapper::SolverParam();
  LPWrapper::SolverParam* snull_ptr = 0;
	TEST_NOT_EQUAL(sptr, snull_ptr)
  TEST_EQUAL(sptr->message_level,3)
  TEST_EQUAL(sptr->branching_tech,4)
  TEST_EQUAL(sptr->backtrack_tech,3)
  TEST_EQUAL(sptr->preprocessing_tech,2)
  TEST_EQUAL(sptr->enable_feas_pump_heuristic,true)
  TEST_EQUAL(sptr->enable_gmi_cuts,true)
  TEST_EQUAL(sptr->enable_mir_cuts,true)
  TEST_EQUAL(sptr->enable_cov_cuts,true)
  TEST_EQUAL(sptr->enable_clq_cuts,true)
  TEST_EQUAL(sptr->mip_gap,0.0)                        
  TEST_EQUAL(sptr->output_freq,5000)
  TEST_EQUAL(sptr->output_delay,10000)
  TEST_EQUAL(sptr->enable_presolve,true)
  TEST_EQUAL(sptr->enable_binarization,true) 
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



