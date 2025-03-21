// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h> // for String

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/config.h>

#include <limits>

// do NOT include glpk and CoinOr headers here, as they define bad stuff, which ripples through OpenMS then...
// include them in LPWrapper.cpp where they do not harm
// only declare them here
class CoinModel;

// if GLPK was found:
#ifndef OPENMS_HAS_COINOR
  #ifndef GLP_PROB_DEFINED
    #define GLP_PROB_DEFINED
    // depending on the glpk version
    // define glp_prob as forward or struct
    #if OPENMS_GLPK_VERSION_MAJOR == 4 && OPENMS_GLPK_VERSION_MINOR < 48
    typedef struct
    {
      double _opaque_prob[100];
    } glp_prob;
    #else
    class glp_prob;
    #endif
  #endif
#endif

namespace OpenMS
{

  class OPENMS_DLLAPI LPWrapper
  {
public:
    /**
       @brief Struct that holds the parameters of the LP solver
    */
    struct SolverParam
    {
      SolverParam() :
        message_level(3), branching_tech(4), backtrack_tech(3),
        preprocessing_tech(2), enable_feas_pump_heuristic(true), enable_gmi_cuts(true),
        enable_mir_cuts(true), enable_cov_cuts(true), enable_clq_cuts(true), mip_gap(0.0),
        time_limit((std::numeric_limits<Int>::max)()), output_freq(5000), output_delay(10000), enable_presolve(true),
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
      double mip_gap;
      Int time_limit;
      Int output_freq;
      Int output_delay;
      bool enable_presolve;
      bool enable_binarization; ///< only with presolve
    };

    enum Type
    {
      UNBOUNDED = 1,
      LOWER_BOUND_ONLY,
      UPPER_BOUND_ONLY,
      DOUBLE_BOUNDED,
      FIXED
    };

    enum VariableType
    {
      CONTINUOUS = 1,
      INTEGER,
      BINARY
    };

    enum Sense
    {
      MIN = 1,
      MAX
    };

    enum WriteFormat
    {
      FORMAT_LP = 0,
      FORMAT_MPS,
      FORMAT_GLPK
    };

    enum SOLVER
    {
      SOLVER_GLPK = 0
#ifdef OPENMS_HAS_COINOR
      , SOLVER_COINOR
#endif
    };

    enum SolverStatus
    {
      UNDEFINED = 1,
      OPTIMAL = 5,
      FEASIBLE = 2,
      NO_FEASIBLE_SOL = 4
    };

    LPWrapper();
    virtual ~LPWrapper();

    // problem creation/manipulation
    /// adds a row to the LP matrix, returns index
    Int addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values, const String& name);
    /// adds an empty column to the LP matrix, returns index
    Int addColumn();
    /// adds a column to the LP matrix, returns index
    Int addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values, const String& name);

    /**
      @brief Adds a row with boundaries to the LP matrix, returns index

      If you have a fixed variable, GLPK requires to use the "fixed" type, instead of "double-bounded" with equal bounds.

      @param row_indices
      @param row_values
      @param name
      @param lower_bound
      @param upper_bound
      @param type Type of the row 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
    */
    Int addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values,
               const String& name, double lower_bound, double upper_bound, Type type);

    /**
      @brief Adds a column with boundaries to the LP matrix, returns index

      @param column_indices
      @param column_values
      @param name
      @param lower_bound
      @param upper_bound
      @param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
    */
    Int addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values, const String& name, double lower_bound, double upper_bound, Type type);

    /// delete index-th row
    void deleteRow(Int index);
    /// sets name of the index-th column
    void setColumnName(Int index, const String& name);
    /// gets name of the index-th column
    String getColumnName(Int index);
    /// sets name of the index-th row
    String getRowName(Int index);
    /// gets index of the row with name
    Int getRowIndex(const String& name);
    /// gets index of the column with name
    Int getColumnIndex(const String& name);
    /// gets column's upper bound
    double getColumnUpperBound(Int index);
    /// gets column's lower bound
    double getColumnLowerBound(Int index);
    /// gets row's upper bound
    double getRowUpperBound(Int index);
    /// gets row's lower bound
    double getRowLowerBound(Int index);
    /// sets name of the index-th row
    void setRowName(Int index, const String& name);

    /**
      @brief Set column bounds.

      @param index
      @param lower_bound
      @param upper_bound
      @param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
     */
    void setColumnBounds(Int index, double lower_bound, double upper_bound, Type type);

    /**
      @brief Set row bounds.

      @param index
      @param lower_bound
      @param upper_bound
      @param type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed constraint
     */
    void setRowBounds(Int index, double lower_bound, double upper_bound, Type type);

    /**
      @brief Set column/variable type.

      @param index
      @param type 1- continuous, 2- integer, 3- binary variable
     */
    void setColumnType(Int index, VariableType type);

    /**
      @brief Get column/variable type.

      @param index
      @return 1- continuous, 2- integer, 3- binary variable
     */
    VariableType getColumnType(Int index);

    /// set objective value for column with index
    void setObjective(Int index, double obj_value);
    /// get objective value for column with index
    double getObjective(Int index);

    /**
      @brief Set objective direction.

      @param sense 1- minimize, 2- maximize
     */
    void setObjectiveSense(Sense sense);
    Sense getObjectiveSense();

    /// get number of columns
    Int getNumberOfColumns();
    /// get number of rows
    Int getNumberOfRows();

    void setElement(Int row_index, Int column_index, double value);
    double getElement(Int row_index, Int column_index);

    // problem reading/writing
    /**
      @brief Read LP from file

      @param filename Filename where to store the LP problem.
      @param format LP, MPS or GLPK.
     */
    void readProblem(const String& filename, const String& format);

    /**
      @brief Write LP formulation to a file.

      @param filename output filename, if the filename ends with '.gz' it will be compressed
      @param format MPS-format is supported by GLPK and COIN-OR; LP and GLPK-formats only by GLPK
     */
    void writeProblem(const String& filename, const WriteFormat format) const;

    /**
      @brief solve problems, parameters like enabled heuristics can be given via solver_param

      The verbose level (0,1,2) determines if the solver prints status messages and internals.

      @param solver_param
      @param verbose_level

      @return solver dependent (todo: fix)
    */
    Int solve(SolverParam& solver_param, const Size verbose_level = 0);

    /**
      @brief Get solution status.

      @return status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution
     */
    SolverStatus getStatus();

    // solution access
    double getObjectiveValue();
    double getColumnValue(Int index);

    Int getNumberOfNonZeroEntriesInRow(Int idx);
    void getMatrixRow(Int idx, std::vector<Int>& indexes);

    /// get currently active solver
    SOLVER getSolver() const;

protected:
#ifdef OPENMS_HAS_COINOR
    CoinModel * model_ = nullptr;
    std::vector<double> solution_;
#else
    glp_prob * lp_problem_ = nullptr;
#endif

    SOLVER solver_;


  }; // class

} // namespace

