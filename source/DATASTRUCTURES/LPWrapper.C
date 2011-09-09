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
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{


  LPWrapper::LPWrapper()
  {
    lp_problem_ = glp_create_prob();
  }

  LPWrapper::~LPWrapper()
  {
  }
  
  Size LPWrapper::addRow(std::vector<Int> row_indices,std::vector<DoubleReal> row_values,String name) // return index
  {
    if(row_indices.size() != row_values.size())  throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"indices and values vectors differ in size");
    Int index = glp_add_rows(lp_problem_, 1);
    // glpk accesses arrays beginning at index 1-> we have to insert an empty value at the front
    row_indices.insert(row_indices.begin(),-1);
    row_values.insert(row_values.begin(),-1);
    for(Size i = 0; i< row_indices.size();++i)   row_indices[i] +=1;//std::cout << row_indices[i]
    glp_set_mat_row(lp_problem_, index, row_indices.size()-1,&(row_indices[0]), &(row_values[0]));
    glp_set_row_name(lp_problem_, index, name.c_str());    
    return index-1;
  }
  
  Size LPWrapper::addColumn()
  {
    return glp_add_cols(lp_problem_, 1)-1;
  }
    
  Size LPWrapper::addColumn(std::vector<Int> column_indices,std::vector<DoubleReal> column_values,String name)
  {
    if(column_indices.size() != column_values.size())  throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"indices and values vectors differ in size");
    Int index = glp_add_cols(lp_problem_, 1);
    // glpk accesses arrays beginning at index 1-> we have to insert an empty value at the front
    column_indices.insert(column_indices.begin(),-1);
    column_values.insert(column_values.begin(),-1);
    for(Size i = 0; i< column_indices.size();++i)   column_indices[i] +=1;
    glp_set_mat_col(lp_problem_, index, column_indices.size()-1,&(column_indices[0]), &(column_values[0]));
    glp_set_col_name(lp_problem_, index, name.c_str());    
    return index-1;
  }

  Size LPWrapper::addRow(std::vector<Int>& row_indices,std::vector<DoubleReal>& row_values,String name,DoubleReal lower_bound,
                         DoubleReal upper_bound,Int type)
  {
    Int index = addRow(row_indices,row_values,name);
    glp_set_row_bnds(lp_problem_,index+1, type, lower_bound, upper_bound);
    return index; // in addRow index is decreased already
  }

  Size LPWrapper::addColumn(std::vector<Int>& column_indices,std::vector<DoubleReal>& column_values,String name,
                            DoubleReal lower_bound,DoubleReal upper_bound,Int type) //return index
  {
    Int index = addColumn(column_indices,column_values,name);
    glp_set_col_bnds(lp_problem_, index+1, type, lower_bound, upper_bound);
    return index;// in addColumn index is decreased already
  }


  void LPWrapper::setColumnName(Size index,String name)
  {
    glp_set_col_name(lp_problem_, (int) index+1, name.c_str());
  }

  void LPWrapper::setRowName(Size index,String name)
  {
    glp_set_row_name(lp_problem_, (int) index+1, name.c_str());
  }

  void LPWrapper::setColumnBounds(Size index,DoubleReal lower_bound,DoubleReal upper_bound,Int type)
  {
    glp_set_col_bnds(lp_problem_, (int) index+1, type, lower_bound, upper_bound);
  }

  void LPWrapper::setRowBounds(Size index,DoubleReal lower_bound,DoubleReal upper_bound,Int type)
  {
    glp_set_row_bnds(lp_problem_, (int) index+1, type, lower_bound, upper_bound);
  }

  void LPWrapper::setColumnType(Size index,Size type) // 0- continuous, 1- integer, 2- binary
  {
    glp_set_col_kind(lp_problem_, (int) index+1, (int) type);
  }

  Int LPWrapper::getColumnType(Size index)
  {
    return glp_get_col_kind(lp_problem_, (int) index+1);
  }

  void LPWrapper::setObjective(Size index,DoubleReal obj_value)
  {
    glp_set_obj_coef(lp_problem_, (int) index+1, obj_value);
  }

  void LPWrapper::setObjectiveSense(Size sense) // 1 min, 2 max
  {
    glp_set_obj_dir(lp_problem_, (int) sense);
  }

  Size LPWrapper::getNumberOfColumns()
  {
    return glp_get_num_cols(lp_problem_);
  }

  Size LPWrapper::getNumberOfRows()
  {
    return glp_get_num_rows(lp_problem_);
  }

  String LPWrapper::getColumnName(Size index)
  {
    return String(glp_get_col_name(lp_problem_, (int) index));
  }
  
  String LPWrapper::getRowName(Size index)
  {
    return String(glp_get_row_name(lp_problem_, (int) index));
  }

  Size LPWrapper::getRowIndex(String name)
  {
    return glp_find_row(lp_problem_, name.c_str());
  }

  Size LPWrapper::getColumnIndex(String name)
  {
    return glp_find_col(lp_problem_, name.c_str());
  }

  void LPWrapper::readProblem(String filename,String format) // format=(LP,MPS,GLPK)
  {
    if(format == "LP")
      {
        glp_read_lp(lp_problem_,NULL,filename.c_str());
      }
    else if(format == "MPS")
      {
        glp_read_mps(lp_problem_,GLP_MPS_FILE,NULL,filename.c_str());
      }
    else if(format == "GLPK")
      {
        glp_read_prob(lp_problem_,0,filename.c_str());
      }
    else  throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"invalid LP format, allowed are LP, MPS, GLPK");
  }

  void LPWrapper::writeProblem(String filename,String format)// format=(LP,MPS,GLPK)
  {
    if(format == "LP")
      {
        glp_write_lp(lp_problem_,NULL,filename.c_str());
      }
    else if(format == "MPS")
      {
        glp_write_mps(lp_problem_,GLP_MPS_FILE,NULL,filename.c_str());
      }
    else if(format == "GLPK")
      {
        glp_write_prob(lp_problem_,0,filename.c_str());
      }
    else  throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"invalid LP format, allowed are LP, MPS, GLPK");

  }


  Int LPWrapper::solve(SolverParam& solver_param) // ruft glp_intopt auf als MIP-Solver, benutzt branch-and-cut
  {
    glp_iocp solver_param_glp;
    glp_init_iocp(&solver_param_glp);		

    solver_param_glp.msg_lev = solver_param.message_level;
    solver_param_glp.br_tech = solver_param.branching_tech;
    solver_param_glp.bt_tech = solver_param.backtrack_tech;
    solver_param_glp.pp_tech = solver_param.preprocessing_tech;
    if(solver_param.enable_feas_pump_heuristic) solver_param_glp.fp_heur = GLP_ON;
    if(solver_param.enable_gmi_cuts) solver_param_glp.gmi_cuts=GLP_ON;
    if(solver_param.enable_mir_cuts) solver_param_glp.mir_cuts=GLP_ON;
    if(solver_param.enable_cov_cuts) solver_param_glp.cov_cuts=GLP_ON;
    if(solver_param.enable_clq_cuts) solver_param_glp.clq_cuts=GLP_ON;
    solver_param_glp.mip_gap = solver_param.mip_gap;
    solver_param_glp.tm_lim = solver_param.time_limit;
    solver_param_glp.out_frq = solver_param.output_freq;
    solver_param_glp.out_dly = solver_param.output_delay;
    if(solver_param.enable_presolve) solver_param_glp.presolve=GLP_ON;
    if(solver_param.enable_binarization) solver_param_glp.binarize=GLP_ON; // only with presolve

    return glp_intopt(lp_problem_, &solver_param_glp);
  }

  Int LPWrapper::getStatus()
  {
    return glp_mip_status(lp_problem_);
  }
  
  DoubleReal LPWrapper::getObjectiveValue()

  {
    return glp_mip_obj_val(lp_problem_);
  }

  DoubleReal LPWrapper::getColumnValue(Size index)
  {
    // glpk uses arrays beginning at pos 1, so we need to shift
    return glp_mip_col_val(lp_problem_, (int) index +1);
  }

  DoubleReal LPWrapper::getRowValue(Size index)
  {
    // glpk uses arrays beginning at pos 1, so we need to shift
    return glp_mip_row_val(lp_problem_, (int) index +1);
  }


}// namespace
