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

// do NOT include glpk, CoinOr, and HIGHS headers here, as they define bad stuff, which ripples through OpenMS then...
// include them in LPWrapper.cpp where they do not harm
// only declare them here
class CoinModel;
class Highs;

// if GLPK was found:
#if !defined(OPENMS_HAS_COINOR) && !defined(OPENMS_HAS_HIGHS)
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

  /**
    @brief A wrapper class for linear programming (LP) solvers

    This class provides a unified interface to different linear programming solvers,
    including GLPK (GNU Linear Programming Kit) and COIN-OR (if available).
    
    Linear programming is a method to find the best outcome in a mathematical model
    whose requirements are represented by linear relationships. It is used for
    optimization problems where the objective function and constraints are linear.
    
    LPWrapper allows you to:
    - Create and manipulate LP problems (add rows, columns, set bounds)
    - Set objective functions and constraints
    - Solve the LP problem using different solvers
    - Access the solution and status information
    
    The class supports both continuous and integer variables, allowing for
    mixed-integer linear programming (MILP) problems.
    
    @ingroup Datastructures
  */
  class OPENMS_DLLAPI LPWrapper
  {
public:
    /**
       @brief Struct that holds the parameters of the LP solver
       
       This structure contains various parameters that control the behavior of the LP solver,
       including algorithm selection, cut generation, heuristics, and output control.
       
       Most parameters have reasonable defaults and don't need to be modified for basic use cases.
       Advanced users can tune these parameters to improve performance for specific problem types.
    */
    struct SolverParam
    {
      /**
        @brief Default constructor that initializes all parameters with reasonable defaults
      */
      SolverParam() :
        message_level(3), branching_tech(4), backtrack_tech(3),
        preprocessing_tech(2), enable_feas_pump_heuristic(true), enable_gmi_cuts(true),
        enable_mir_cuts(true), enable_cov_cuts(true), enable_clq_cuts(true), mip_gap(0.0),
        time_limit((std::numeric_limits<Int>::max)()), output_freq(5000), output_delay(10000), enable_presolve(true),
        enable_binarization(true)
      {
      }

      Int message_level;                ///< Controls verbosity of solver output (0-3)
      Int branching_tech;               ///< Branching technique for MIP problems
      Int backtrack_tech;               ///< Backtracking technique for MIP problems
      Int preprocessing_tech;           ///< Preprocessing technique
      bool enable_feas_pump_heuristic;  ///< Enable feasibility pump heuristic
      bool enable_gmi_cuts;             ///< Enable Gomory mixed-integer cuts
      bool enable_mir_cuts;             ///< Enable mixed-integer rounding cuts
      bool enable_cov_cuts;             ///< Enable cover cuts
      bool enable_clq_cuts;             ///< Enable clique cuts
      double mip_gap;                   ///< Relative gap tolerance for MIP problems
      Int time_limit;                   ///< Time limit in milliseconds
      Int output_freq;                  ///< Output frequency in milliseconds
      Int output_delay;                 ///< Output delay in milliseconds
      bool enable_presolve;             ///< Enable presolve techniques
      bool enable_binarization;         ///< Enable binarization (only with presolve)
    };

    /**
      @brief Enumeration for variable/constraint bound types
      
      Defines the type of bounds applied to variables or constraints in the LP problem.
    */
    enum Type
    {
      UNBOUNDED = 1,      ///< No bounds (free variable)
      LOWER_BOUND_ONLY,   ///< Only lower bound is specified
      UPPER_BOUND_ONLY,   ///< Only upper bound is specified
      DOUBLE_BOUNDED,     ///< Both lower and upper bounds are specified
      FIXED               ///< Lower bound equals upper bound (fixed value)
    };

    /**
      @brief Enumeration for variable types in the LP problem
      
      Defines whether variables are continuous or discrete (integer/binary).
    */
    enum VariableType
    {
      CONTINUOUS = 1,     ///< Continuous variable (can take any real value within bounds)
      INTEGER,            ///< Integer variable (can only take integer values within bounds)
      BINARY              ///< Binary variable (can only take values 0 or 1)
    };

    /**
      @brief Enumeration for optimization direction
      
      Defines whether the objective function should be minimized or maximized.
    */
    enum Sense
    {
      MIN = 1,            ///< Minimize the objective function
      MAX                 ///< Maximize the objective function
    };

    /**
      @brief Enumeration for LP problem file formats
      
      Defines the file format used when writing LP problems to disk.
    */
    enum WriteFormat
    {
      FORMAT_LP = 0,      ///< LP format (human-readable)
      FORMAT_MPS,         ///< MPS format (industry standard)
      FORMAT_GLPK         ///< GLPK's native format
    };

    /**
      @brief Enumeration for available LP solvers
      
      Defines which solver backend to use for solving LP problems.
    */
    enum SOLVER
    {
      SOLVER_GLPK = 0     ///< GNU Linear Programming Kit solver
#ifdef OPENMS_HAS_COINOR
      , SOLVER_COINOR     ///< COIN-OR solver (if available)
#endif
#ifdef OPENMS_HAS_HIGHS
      , SOLVER_HIGHS      ///< HIGHS solver (if available)
#endif
    };

    /**
      @brief Enumeration for solver status after solving an LP problem
      
      Defines the possible outcomes after attempting to solve an LP problem.
    */
    enum SolverStatus
    {
      UNDEFINED = 1,      ///< Status is undefined (e.g., solver not run yet)
      OPTIMAL = 5,        ///< Optimal solution found
      FEASIBLE = 2,       ///< Feasible solution found (but not necessarily optimal)
      NO_FEASIBLE_SOL = 4 ///< No feasible solution exists for the problem
    };

    /**
      @brief Default constructor
      
      Initializes a new LP problem with the default solver (GLPK or COIN-OR if available).
    */
    LPWrapper();
    
    /**
      @brief Virtual destructor
      
      Frees all resources associated with the LP problem.
    */
    virtual ~LPWrapper();

    // problem creation/manipulation
    /**
      @brief Adds a row to the LP matrix
      
      @param[in] row_indices Indices of the columns that have non-zero coefficients in this row
      @param[in] row_values Values of the non-zero coefficients in this row
      @param[in] name Name of the row (for identification purposes)
      @return Index of the newly added row
    */
    Int addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values, const String& name);
    
    /**
      @brief Adds an empty column to the LP matrix
      
      @return Index of the newly added column
    */
    Int addColumn();
    
    /**
      @brief Adds a column to the LP matrix
      
      @param[in] column_indices Indices of the rows that have non-zero coefficients in this column
      @param[in] column_values Values of the non-zero coefficients in this column
      @param[in] name Name of the column (for identification purposes)
      @return Index of the newly added column
    */
    Int addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values, const String& name);

    /**
      @brief Adds a row with boundaries to the LP matrix, returns index

      If you have a fixed variable, GLPK requires to use the "fixed" type, instead of "double-bounded" with equal bounds.

      @param[in] row_indices
      @param[in] row_values
      @param[in] name
      @param[in] lower_bound
      @param[in] upper_bound
      @param[in] type Type of the row 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
    */
    Int addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values,
               const String& name, double lower_bound, double upper_bound, Type type);

    /**
      @brief Adds a column with boundaries to the LP matrix, returns index

      @param[in] column_indices
      @param[in] column_values
      @param[in] name
      @param[in] lower_bound
      @param[in] upper_bound
      @param[in] type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
    */
    Int addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values, const String& name, double lower_bound, double upper_bound, Type type);

    /**
      @brief Delete the row at the specified index
      
      @param[in] index Index of the row to delete
    */
    void deleteRow(Int index);
    
    /**
      @brief Set the name of a column
      
      @param[in] index Index of the column to rename
      @param[in] name New name for the column
    */
    void setColumnName(Int index, const String& name);
    
    /**
      @brief Get the name of a column
      
      @param[in] index Index of the column
      @return Name of the column
    */
    String getColumnName(Int index);
    
    /**
      @brief Get the name of a row
      
      @param[in] index Index of the row
      @return Name of the row
    */
    String getRowName(Int index);
    
    /**
      @brief Find the index of a row by its name
      
      @param[in] name Name of the row to find
      @return Index of the row with the given name
    */
    Int getRowIndex(const String& name);
    
    /**
      @brief Find the index of a column by its name
      
      @param[in] name Name of the column to find
      @return Index of the column with the given name
    */
    Int getColumnIndex(const String& name);
    
    /**
      @brief Get the upper bound of a column
      
      @param[in] index Index of the column
      @return Upper bound value of the column
    */
    double getColumnUpperBound(Int index);
    
    /**
      @brief Get the lower bound of a column
      
      @param[in] index Index of the column
      @return Lower bound value of the column
    */
    double getColumnLowerBound(Int index);
    
    /**
      @brief Get the upper bound of a row
      
      @param[in] index Index of the row
      @return Upper bound value of the row
    */
    double getRowUpperBound(Int index);
    
    /**
      @brief Get the lower bound of a row
      
      @param[in] index Index of the row
      @return Lower bound value of the row
    */
    double getRowLowerBound(Int index);
    
    /**
      @brief Set the name of a row
      
      @param[in] index Index of the row to rename
      @param[in] name New name for the row
    */
    void setRowName(Int index, const String& name);

    /**
      @brief Set column bounds.

      @param[in] index
      @param[in] lower_bound
      @param[in] upper_bound
      @param[in] type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed variable
     */
    void setColumnBounds(Int index, double lower_bound, double upper_bound, Type type);

    /**
      @brief Set row bounds.

      @param[in] index
      @param[in] lower_bound
      @param[in] upper_bound
      @param[in] type 1 - unbounded, 2 - only lower bound, 3 - only upper bound, 4 - double-bounded variable, 5 - fixed constraint
     */
    void setRowBounds(Int index, double lower_bound, double upper_bound, Type type);

    /**
      @brief Set column/variable type.

      @param[in] index
      @param[in] type 1- continuous, 2- integer, 3- binary variable
     */
    void setColumnType(Int index, VariableType type);

    /**
      @brief Get column/variable type.

      @param[in] index
      @return 1- continuous, 2- integer, 3- binary variable
     */
    VariableType getColumnType(Int index);

    /**
      @brief Set the objective coefficient for a column/variable
      
      @param[in] index Index of the column/variable
      @param[in] obj_value Coefficient value in the objective function
    */
    void setObjective(Int index, double obj_value);
    
    /**
      @brief Get the objective coefficient for a column/variable
      
      @param[in] index Index of the column/variable
      @return Coefficient value in the objective function
    */
    double getObjective(Int index);

    /**
      @brief Set objective direction.

      @param[in] sense 1- minimize, 2- maximize
     */
    void setObjectiveSense(Sense sense);
    /**
      @brief Get the current objective direction
      
      @return Current optimization direction (MIN or MAX)
    */
    Sense getObjectiveSense();

    /**
      @brief Get the number of columns/variables in the LP problem
      
      @return Number of columns in the LP matrix
    */
    Int getNumberOfColumns();
    
    /**
      @brief Get the number of rows/constraints in the LP problem
      
      @return Number of rows in the LP matrix
    */
    Int getNumberOfRows();

    /**
      @brief Set the value of a matrix element at the specified position
      
      @param[in] row_index Index of the row
      @param[in] column_index Index of the column
      @param[in] value Value to set at the specified position
    */
    void setElement(Int row_index, Int column_index, double value);
    
    /**
      @brief Get the value of a matrix element at the specified position
      
      @param[in] row_index Index of the row
      @param[in] column_index Index of the column
      @return Value at the specified position
    */
    double getElement(Int row_index, Int column_index);

    // problem reading/writing
    /**
      @brief Read LP from file

      @param[out] filename Filename where to store the LP problem.
      @param[in] format LP, MPS or GLPK.
     */
    void readProblem(const String& filename, const String& format);

    /**
      @brief Write LP formulation to a file.

      @param[out] filename output filename, if the filename ends with '.gz' it will be compressed
      @param[in] format MPS-format is supported by GLPK and COIN-OR; LP and GLPK-formats only by GLPK
     */
    void writeProblem(const String& filename, const WriteFormat format) const;

    /**
      @brief solve problems, parameters like enabled heuristics can be given via solver_param

      The verbose level (0,1,2) determines if the solver prints status messages and internals.

      @param[in] solver_param
      @param[in] verbose_level

      @return solver dependent (todo: fix)
    */
    Int solve(SolverParam& solver_param, const Size verbose_level = 0);

    /**
      @brief Get solution status.

      @return status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution
     */
    SolverStatus getStatus();

    // solution access
    /**
      @brief Get the objective function value of the solution
      
      @return Value of the objective function at the optimal solution
    */
    double getObjectiveValue();
    
    /**
      @brief Get the value of a variable in the solution
      
      @param[in] index Index of the column/variable
      @return Value of the variable in the optimal solution
    */
    double getColumnValue(Int index);

    /**
      @brief Get the number of non-zero entries in a specific row
      
      @param[in] idx Index of the row
      @return Number of non-zero coefficients in the row
    */
    Int getNumberOfNonZeroEntriesInRow(Int idx);
    
    /**
      @brief Get the indices of non-zero entries in a specific row
      
      @param[in] idx Index of the row
      @param[out] indexes Vector to store the column indices of non-zero entries
    */
    void getMatrixRow(Int idx, std::vector<Int>& indexes);

    /**
      @brief Get the currently active solver backend
      
      @return Currently active solver (GLPK or COIN-OR)
    */
    SOLVER getSolver() const;

protected:
#ifdef OPENMS_HAS_COINOR
    CoinModel * model_ = nullptr;      ///< COIN-OR model object for the LP problem
    std::vector<double> solution_;     ///< Solution vector when using COIN-OR
#elif defined(OPENMS_HAS_HIGHS)
    Highs * highs_model_ = nullptr;    ///< HIGHS model object for the LP problem
#else
    glp_prob * lp_problem_ = nullptr;  ///< GLPK problem object for the LP problem
#endif

    SOLVER solver_;                    ///< Currently active solver backend


  }; // class

} // namespace

