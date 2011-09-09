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
#ifndef OPENMS_DATASTRUCTURES_LPWRAPPER_H
#define OPENMS_DATASTRUCTURES_LPWRAPPER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <glpk.h>
//#include <limits>
namespace OpenMS
{

  class OPENMS_DLLAPI LPWrapper
  {
  public:
    struct SolverParam
    {
      SolverParam(): message_level(GLP_MSG_ALL),branching_tech(GLP_BR_DTH),backtrack_tech(GLP_BT_BLB),
                     preprocessing_tech(GLP_PP_ALL),enable_feas_pump_heuristic(true),enable_gmi_cuts(true),
                     enable_mir_cuts(true),enable_cov_cuts(true),enable_clq_cuts(true),mip_gap(0.0),
                     time_limit(std::numeric_limits<Int>::max()), output_freq(5000),output_delay(10000),enable_presolve(true),
                     enable_binarization(true)

      {

      }
      Int message_level;
      Int branching_tech;
      Int backtrack_tech;
      Int preprocessing_tech;
      bool enable_feas_pump_heuristic;
      bool enable_gmi_cuts;
      bool enable_mir_cuts;
      bool enable_cov_cuts;
      bool enable_clq_cuts;
      DoubleReal mip_gap;
      Int time_limit;
      Int output_freq;
      Int output_delay;
      bool enable_presolve;
      bool enable_binarization; // only with presolve
    };


    
		LPWrapper();
    virtual ~LPWrapper();

    // problem creation/manipulation
    /// adds a row to the LP matrix, returns index
    Size addRow(std::vector<Int> row_indices,std::vector<DoubleReal> row_values,String name); 
    /// adds an empty column to the LP matrix, returns index
    Size addColumn();
    /// adds a column to the LP matrix, returns index
    Size addColumn(std::vector<Int> column_indices,std::vector<DoubleReal> column_values,String name);
    /*
      @brief Adds a row with boundaries to the LP matrix, returns index
      @param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
    */
    Size addRow(std::vector<Int>& row_indices,std::vector<DoubleReal>& row_values,String name,DoubleReal lower_bound,DoubleReal upper_bound,Int type);
    /*
      @brief Adds a column with boundaries to the LP matrix, returns index
      @param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
    */
    Size addColumn(std::vector<Int>& column_indices,std::vector<DoubleReal>& column_values,String name,DoubleReal lower_bound,DoubleReal upper_bound,Int type);
    /// sets name of the index-th column
    void setColumnName(Size index,String name);
    /// gets name of the index-th column
    String getColumnName(Size index);
    /// sets name of the index-th row
    String getRowName(Size index);
    /// gets index of the row with name
    Size getRowIndex(String name);
    /// gets index of the column with name
    Size getColumnIndex(String name);
    /// sets name of the index-th row
    void setRowName(Size index,String name);
    /**
     *	@brief Set column bounds.
     *	
     *	@param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
     */
    void setColumnBounds(Size index,DoubleReal lower_bound,DoubleReal upper_bound,Int type);
    /**
     *	@brief Set row bounds.
     *	
     *	@param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed constraint
     */
    void setRowBounds(Size index,DoubleReal lower_bound,DoubleReal upper_bound,Int type);
    /**
     *	@brief Set column/variable type.
     *	
     *	@param type 1- continuous, 2- integer, 3- binary variable
     */
    void setColumnType(Size index,Size type);
    /**
     *	@brief Get column/variable type.
     *	
     *	@return 1- continuous, 2- integer, 3- binary variable
     */
    Int getColumnType(Size index);
    /// set objective value for column with index
    void setObjective(Size index,DoubleReal obj_value);
    /**
     *	@brief Set objective direction.
     *	
     *	@param sense 1- minimize, 2- maximize
     */
    void setObjectiveSense(Size sense);
    /// get number of columns
    Size getNumberOfColumns();
    /// get number of rows
    Size getNumberOfRows();
    
    // problem reading/writing
    /**
     *	@brief Read LP from file
     *	
     *	@param format LP, MPS or GLPK
     */
    void readProblem(String filename,String format); 
    /**
     *	@brief Write LP formulation to a file.
     *	
     *	@param filename output filename, if the filename ends with '.gz' it will be compressed
		 *  @param format can be LP, MPS or GLPK
     */
    void writeProblem(String filename,String format);// format=(LP,MPS,GLPK)

    // problem solving
    /// solve problems, parameters like enabled heuristics can be given via solver_param
    Int solve(SolverParam& solver_param); 
    Int getStatus();
    // solution access
    DoubleReal getObjectiveValue();
    DoubleReal getColumnValue(Size index);
    DoubleReal getRowValue(Size index);
	protected:
    glp_prob* lp_problem_;

  }; // class



} // namespace

#endif // OPENMS_DATASTRUCTURES_LPWRAPPER_H
