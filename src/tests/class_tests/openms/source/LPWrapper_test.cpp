// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#ifdef OPENMS_HAS_COINOR
  #ifdef OPENMS_HAS_COIN_INCLUDE_SUBDIR_IS_COIN
    #include "coin/CoinModel.hpp"
  #else
    #include "coin-or/CoinModel.hpp"
  #endif
#endif
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LPWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LPWrapper* ptr = nullptr;
LPWrapper* null_ptr = nullptr;
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
std::vector<double> values(2,0.5);
std::vector<Int> indices;
indices.push_back(0);
indices.push_back(1);
START_SECTION((Int addColumn()))
{
  lp.addColumn();
  TEST_EQUAL(lp.getNumberOfColumns(),1);  
}
END_SECTION



START_SECTION((Int addRow(fewer args)))
{
  lp.addColumn();
  lp.addRow(indices,values,String("row1"));
  TEST_EQUAL(lp.getNumberOfRows(),1);
  TEST_EQUAL(lp.getRowName(0),"row1");
}
END_SECTION


START_SECTION((Int addColumn(std::vector< Int > column_indices, std::vector< double > column_values, const String &name)))
{
  lp.addRow(indices,values,String("row2"));
  lp.addColumn(indices,values,String("col3"));
  TEST_EQUAL(lp.getNumberOfColumns(),3);
  TEST_EQUAL(lp.getColumnName(2),"col3");
}
END_SECTION

START_SECTION((Int addRow(all args)))
{
  lp.addRow(indices,values,String("row3"),0.2,1.2,LPWrapper::DOUBLE_BOUNDED);
  TEST_EQUAL(lp.getNumberOfRows(),3);
  TEST_EQUAL(lp.getRowName(2),"row3");
}
END_SECTION

START_SECTION((Int addColumn(std::vector< Int > &column_indices, std::vector< double > &column_values, const String &name, double lower_bound, double upper_bound, Type type)))
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

START_SECTION((void setColumnBounds(Int index, double lower_bound, double upper_bound, Type type)))
{
  lp.setColumnBounds(0,0.3,1.0,LPWrapper::DOUBLE_BOUNDED);
  TEST_EQUAL(lp.getColumnUpperBound(0),1.0);
  TEST_EQUAL(lp.getColumnLowerBound(0),0.3);
}
END_SECTION

START_SECTION((void setRowBounds(Int index, double lower_bound, double upper_bound, Type type)))
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

START_SECTION((void setObjective(Int index, double obj_value)))
{
  lp.setObjective(0,3.5);
  TEST_EQUAL(lp.getObjective(0),3.5);
}
END_SECTION

START_SECTION((double getObjective(Int index)))
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

START_SECTION((double getColumnUpperBound(Int index)))
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
#ifdef OPENMS_HAS_COINOR
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

START_SECTION((double getColumnLowerBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getColumnLowerBound(0),0.3)
}
END_SECTION

START_SECTION((double getRowUpperBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getRowUpperBound(0),1.0)
}
END_SECTION

START_SECTION((double getRowLowerBound(Int index)))
{
  TEST_REAL_SIMILAR(lp.getRowLowerBound(0),-0.3)  
}
END_SECTION


START_SECTION((void setElement(Int row_index, Int column_index, double value)))
{
  lp.setElement(1,2,0.5);
  TEST_REAL_SIMILAR(lp.getElement(1,2),0.5)  
}
END_SECTION


START_SECTION((double getElement(Int row_index, Int column_index)))
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
#ifdef OPENMS_HAS_COINOR
  else  if (lp.getSolver()==LPWrapper::SOLVER_COINOR)
    {
      lp.readProblem(OPENMS_GET_TEST_DATA_PATH("LPWrapper_test.mps"),"MPS");
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
#endif
}
END_SECTION

START_SECTION((void writeProblem(const String &filename, const WriteFormat format) const ))
{
#ifdef OPENMS_HAS_COINOR
    String tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    lp.writeProblem(tmp_filename,LPWrapper::FORMAT_MPS);
    LPWrapper lp2;
    lp2.readProblem(tmp_filename,"MPS");
    TEST_EQUAL(lp2.getNumberOfColumns(),2)
    TEST_EQUAL(lp2.getNumberOfRows(),3)
    TEST_EQUAL(lp2.getColumnType(0),LPWrapper::INTEGER)
    TEST_EQUAL(lp2.getColumnType(1),LPWrapper::INTEGER)
    TEST_EQUAL(lp2.getObjective(0),1)
    TEST_EQUAL(lp2.getObjective(1),0)
    TEST_EQUAL(lp2.getRowUpperBound(0),0)
    TEST_EQUAL(lp2.getRowUpperBound(1),12)
    TEST_EQUAL(lp2.getRowUpperBound(2),12)
    TEST_EQUAL(lp2.getElement(0,0),1)
    TEST_EQUAL(lp2.getElement(0,1),-1)
    TEST_EQUAL(lp2.getElement(1,0),2)
    TEST_EQUAL(lp2.getElement(1,1),3)
    TEST_EQUAL(lp2.getElement(2,0),3)
    TEST_EQUAL(lp2.getElement(2,1),2)
#else
    String tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    lp.writeProblem(tmp_filename, LPWrapper::FORMAT_LP);
    LPWrapper lp2;
    lp2.readProblem(tmp_filename, "LP");
    TEST_EQUAL(lp2.getNumberOfColumns(), 2)
    TEST_EQUAL(lp2.getNumberOfRows(), 3)
    TEST_EQUAL(lp2.getColumnType(0), LPWrapper::INTEGER)
    TEST_EQUAL(lp2.getColumnType(1), LPWrapper::INTEGER)
    TEST_EQUAL(lp2.getObjective(0), 1)
    TEST_EQUAL(lp2.getObjective(1), 0)
    TEST_EQUAL(lp2.getRowUpperBound(0), 0)
    TEST_EQUAL(lp2.getRowUpperBound(1), 12)
    TEST_EQUAL(lp2.getRowUpperBound(2), 12)
    TEST_EQUAL(lp2.getElement(0, 0), 1)
    TEST_EQUAL(lp2.getElement(0, 1), -1)
    TEST_EQUAL(lp2.getElement(1, 0), 2)
    TEST_EQUAL(lp2.getElement(1, 1), 3)
    TEST_EQUAL(lp2.getElement(2, 0), 3)
    TEST_EQUAL(lp2.getElement(2, 1), 2)
#endif
}
END_SECTION

START_SECTION((Int solve(SolverParam &solver_param, const Size verbose_level=0)))
{
  LPWrapper lp2;
  lp2.readProblem(OPENMS_GET_TEST_DATA_PATH("LPWrapper_test.mps"),"MPS");
  lp2.setObjectiveSense(LPWrapper::MAX);
  LPWrapper::SolverParam param2;
  lp2.solve(param2);
  TEST_EQUAL(lp2.getColumnValue(0),1)
  TEST_EQUAL(lp2.getColumnValue(1),1)

  // Test an integer problem
  LPWrapper  lp3;
  lp3.readProblem(OPENMS_GET_TEST_DATA_PATH("LPWrapper_test_integer.mps"),"MPS");
  lp3.setObjectiveSense(LPWrapper::MAX);
  LPWrapper::SolverParam param3;
  lp3.solve(param3);
  TEST_EQUAL(lp3.getColumnValue(0),2)
  TEST_EQUAL(lp3.getColumnValue(1),2)
}
END_SECTION



// Test an integer problem
LPWrapper  lp4;
lp4.readProblem(OPENMS_GET_TEST_DATA_PATH("LPWrapper_test_integer.mps"),"MPS");
lp4.setObjectiveSense(LPWrapper::MAX);
LPWrapper::SolverParam param4;
lp4.solve(param4);
START_SECTION((SolverStatus getStatus()))
{
  if(lp4.getSolver() == LPWrapper::SOLVER_GLPK)
    {
      TEST_EQUAL(lp4.getStatus(),LPWrapper::OPTIMAL)
    }
#ifdef OPENMS_HAS_COINOR
  else
  {
    TEST_EQUAL(lp4.getStatus(),LPWrapper::UNDEFINED)
  }
#endif

}
END_SECTION

START_SECTION((double getObjectiveValue()))
{
  TEST_REAL_SIMILAR(lp4.getObjectiveValue(),2)
}
END_SECTION

START_SECTION((double getColumnValue(Int index)))
{
  TEST_REAL_SIMILAR(lp4.getColumnValue(0),2)
  TEST_REAL_SIMILAR(lp4.getColumnValue(1),2)
}
END_SECTION

START_SECTION((Int getNumberOfNonZeroEntriesInRow(Int idx)))
{
  TEST_EQUAL(lp4.getNumberOfNonZeroEntriesInRow(0),2)
}
END_SECTION

START_SECTION((void getMatrixRow(Int idx,std::vector<Int>& indexes)))
{
  std::vector<Int> idxs,idxs2;
  idxs.push_back(0);
  idxs.push_back(1);
  lp4.getMatrixRow(0,idxs2);
  TEST_EQUAL(idxs2.size(),idxs.size())
  for(Size i = 0; i < idxs2.size();++i)
    {
      TEST_EQUAL(idxs2[i],idxs[i])
    }
}
END_SECTION

START_SECTION((SOLVER getSolver() const ))
{

#ifdef OPENMS_HAS_COINOR
  TEST_EQUAL(lp4.getSolver(),LPWrapper::SOLVER_COINOR)
#else
  TEST_EQUAL(lp4.getSolver(), LPWrapper::SOLVER_GLPK)
#endif
}
END_SECTION

START_SECTION(([LPWrapper::SolverParam] SolverParam()))
{
  LPWrapper::SolverParam* sptr = new LPWrapper::SolverParam();
  LPWrapper::SolverParam* snull_ptr = nullptr;
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
  delete sptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



