// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>

#ifdef _MSC_VER //disable some COIN-OR warnings that distract from ours
#pragma warning( push ) // save warning state
#pragma warning( disable : 4267 )
#else
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "coin/CoinModel.hpp"
#include "coin/OsiClpSolverInterface.hpp"
#include "coin/CbcModel.hpp"
#include "coin/CbcHeuristic.hpp"
#include "coin/CbcHeuristicLocal.hpp"
#include "coin/CglGomory.hpp"
#include "coin/CglKnapsackCover.hpp"
#include "coin/CglOddHole.hpp"
#include "coin/CglClique.hpp"
#include "coin/CglMixedIntegerRounding.hpp"
#ifdef _MSC_VER
#pragma warning( pop ) // restore old warning state
#else
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif



#include <glpk.h>

namespace OpenMS
{


  LPWrapper::LPWrapper()
  {
#if COINOR_SOLVER == 1
    solver_ = SOLVER_COINOR;
    model_ = new CoinModel;
#else
    solver_ = SOLVER_GLPK;
    lp_problem_ = glp_create_prob();
#endif
  }

  LPWrapper::~LPWrapper()
  {
#if COINOR_SOLVER == 1
    delete model_;
#else
    glp_delete_prob(lp_problem_);
#endif
  }

  Int LPWrapper::addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values, const String& name) // return index
  {
    if (row_indices.size() != row_values.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Indices and values vectors differ in size");
    }
#if COINOR_SOLVER == 1
    model_->addRow(row_indices.size(), row_indices.data(), row_values.data(), -COIN_DBL_MAX, COIN_DBL_MAX, name.c_str());
    return model_->numberRows() - 1;
#else
    std::vector<Int> row_indices_ = row_indices;
    std::vector<double> row_values_ = row_values;
    const Int index = glp_add_rows(lp_problem_, 1);
    // glpk accesses arrays beginning at index 1-> we have to insert an empty value at the front
    row_indices_.insert(row_indices_.begin(), -1);
    row_values_.insert(row_values_.begin(), -1);
    for (Int& row_index : row_indices_)
    {
      ++row_index;
    }
    glp_set_mat_row(lp_problem_, index, row_indices_.size() - 1, row_indices_.data(), row_values_.data());
    glp_set_row_name(lp_problem_, index, name.c_str());
    return index - 1;
#endif
  }

  Int LPWrapper::addColumn()
  {
#if COINOR_SOLVER == 1
    model_->addColumn(0, nullptr, nullptr, 0, 0); // new columns are initially fixed at zero, like in glpk
    return model_->numberColumns() - 1;
#else
    return glp_add_cols(lp_problem_, 1) - 1;
#endif
  }

  Int LPWrapper::addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values, const String& name)
  {
    if (column_indices.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Column indices for Row are empty");
    }
    if (column_indices.size() != column_values.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Indices and values vectors differ in size");
    }
#if COINOR_SOLVER == 1
    model_->addColumn(column_indices.size(), column_indices.data(), column_values.data(), -COIN_DBL_MAX, COIN_DBL_MAX, 0.0, name.c_str());
    return model_->numberColumns() - 1;
#else
    std::vector<Int> column_indices_ = column_indices;
    std::vector<double> column_values_ = column_values;
    const Int index = glp_add_cols(lp_problem_, 1);
    // glpk accesses arrays beginning at index 1-> we have to insert an empty value at the front
    column_indices_.insert(column_indices_.begin(), -1);
    column_values_.insert(column_values_.begin(), -1);
    for (Int& column_index : column_indices_)
    {
      ++column_index;
    }
    glp_set_mat_col(lp_problem_, index, column_indices_.size() - 1, column_indices_.data(), column_values_.data());
    glp_set_col_name(lp_problem_, index, name.c_str());
    return index - 1;
#endif
  }

  Int LPWrapper::addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values,
                        const String& name, double lower_bound, double upper_bound, Type type)
  {
    const Int index = addRow(row_indices, row_values, name);
#if COINOR_SOLVER == 1
    switch (type)
    {
    case UNBOUNDED: // unbounded
      model_->setRowBounds(index, -COIN_DBL_MAX, COIN_DBL_MAX);
      break;

    case LOWER_BOUND_ONLY: // only lower bound
      model_->setRowBounds(index, lower_bound, COIN_DBL_MAX);
      break;

    case UPPER_BOUND_ONLY: // only upper bound
      model_->setRowBounds(index, -COIN_DBL_MAX, upper_bound);
      break;

    default: // double-bounded or fixed
      model_->setRowBounds(index, lower_bound, upper_bound);
      break;
    }
#else
    glp_set_row_bnds(lp_problem_, index + 1, type, lower_bound, upper_bound);
#endif
    return index; // in addRow index is decreased already
  }

  Int LPWrapper::addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values,
                           const String& name, double lower_bound, double upper_bound, Type type) //return index
  {
    const Int index = addColumn(column_indices, column_values, name);
#if COINOR_SOLVER == 1
    switch (type)
    {
    case UNBOUNDED: // unbounded
      model_->setColumnBounds(index, -COIN_DBL_MAX, COIN_DBL_MAX);
      break;

    case LOWER_BOUND_ONLY: // only lower bound
      model_->setColumnBounds(index, lower_bound, COIN_DBL_MAX);
      break;

    case UPPER_BOUND_ONLY: // only upper bound
      model_->setColumnBounds(index, -COIN_DBL_MAX, upper_bound);
      break;

    default: // double-bounded or fixed
      model_->setColumnBounds(index, lower_bound, upper_bound);
      break;
    }
#else
    glp_set_col_bnds(lp_problem_, index + 1, type, lower_bound, upper_bound);
#endif
    return index; // in addColumn index is decreased already
  }

  void LPWrapper::deleteRow(Int index)
  {
#if COINOR_SOLVER == 1
    model_->deleteRow(index);
#else
    int num[] = { 0, index + 1 }; // glpk starts reading at pos 1
    glp_del_rows(lp_problem_, 1, num);
#endif
  }

  void LPWrapper::setElement(Int row_index, Int column_index, double value)
  {
    if (row_index >= getNumberOfRows() || column_index >= getNumberOfColumns())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid index given", String("invalid column_index or row_index"));
    }
#if COINOR_SOLVER == 1
    model_->setElement(row_index, column_index, value);
#else
    const Int length = glp_get_mat_row(lp_problem_, row_index + 1, nullptr, nullptr); // get row length
    std::vector<double> values(length + 1);
    std::vector<Int> indices(length + 1);
    glp_get_mat_row(lp_problem_, row_index + 1, indices.data(), values.data());
    bool found = false;
    for (Int i = 1; i <= length; ++i)
    {
      if (indices[i] == column_index + 1)
      {
        values[i] = value;
        found = true;
        break;
      }
    }
    if (!found) // if this entry wasn't existing before we have to enter it
    {
      std::vector<Int> n_indices(length + 2);
      std::vector<double> n_values(length + 2);
      for (Int i = 0; i <= length; ++i)
      {
        n_indices[i] = indices[i];
        n_values[i] = values[i];
      }
      // now add new value
      n_indices[length + 1] = column_index + 1; // glpk starts reading at pos 1
      n_values[length + 1] = value;
      glp_set_mat_row(lp_problem_, row_index + 1, length, n_indices.data(), n_values.data());
    }
    else
    {
      glp_set_mat_row(lp_problem_, row_index + 1, length, indices.data(), values.data());
    }
#endif
  }

  double LPWrapper::getElement(Int row_index, Int column_index)
  {
    if (row_index >= getNumberOfRows() || column_index >= getNumberOfColumns())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid index given", String("invalid column_index or row_index"));
    }
#if COINOR_SOLVER == 1
    return model_->getElement(row_index, column_index);
#else
    const Int length = glp_get_mat_row(lp_problem_, row_index + 1, nullptr, nullptr);
    std::vector<double> values(length + 1);
    std::vector<Int> indices(length + 1);
    glp_get_mat_row(lp_problem_, row_index + 1, indices.data(), values.data());
    for (Int i = 1; i <= length; ++i)
    {
      if (indices[i] == column_index + 1)
        return values[i];
    }
    return 0.;
#endif
  }

  void LPWrapper::setColumnName(Int index, const String& name)
  {
#if COINOR_SOLVER == 1
    model_->setColumnName(index, name.c_str());
#else
    glp_set_col_name(lp_problem_, index + 1, name.c_str());
#endif
  }

  void LPWrapper::setRowName(Int index, const String& name)
  {
#if COINOR_SOLVER == 1
    model_->setRowName(index, name.c_str());
#else
    glp_set_row_name(lp_problem_, index + 1, name.c_str());
#endif
  }

  void LPWrapper::setColumnBounds(Int index, double lower_bound, double upper_bound, Type type)
  {
#if COINOR_SOLVER == 1
    switch (type)
    {
    case UNBOUNDED: // unbounded
      model_->setColumnBounds(index, -COIN_DBL_MAX, COIN_DBL_MAX);
      break;

    case LOWER_BOUND_ONLY: // only lower bound
      model_->setColumnBounds(index, lower_bound, COIN_DBL_MAX);
      break;

    case UPPER_BOUND_ONLY: // only upper bound
      model_->setColumnBounds(index, -COIN_DBL_MAX, upper_bound);
      break;

    default: // double-bounded or fixed
      model_->setColumnBounds(index, lower_bound, upper_bound);
      break;
    }
#else
    glp_set_col_bnds(lp_problem_, index + 1, type, lower_bound, upper_bound);
#endif
  }

  void LPWrapper::setRowBounds(Int index, double lower_bound, double upper_bound, Type type)
  {
#if COINOR_SOLVER == 1
    switch (type)
    {
    case UNBOUNDED: // unbounded
      model_->setRowBounds(index, -COIN_DBL_MAX, COIN_DBL_MAX);
      break;

    case LOWER_BOUND_ONLY: // only lower bound
      model_->setRowBounds(index, lower_bound, COIN_DBL_MAX);
      break;

    case UPPER_BOUND_ONLY: // only upper bound
      model_->setRowBounds(index, -COIN_DBL_MAX, upper_bound);
      break;

    default: // double-bounded or fixed
      model_->setRowBounds(index, lower_bound, upper_bound);
      break;
    }
#else
    glp_set_row_bnds(lp_problem_, index + 1, type, lower_bound, upper_bound);
#endif
  }

  void LPWrapper::setColumnType(Int index, VariableType type) // 1- continuous, 2- integer, 3- binary
  {
#if COINOR_SOLVER == 1
    if (type == 1)
      model_->setContinuous(index);
    else if (type == 3)
    {
      OPENMS_LOG_WARN << "Coin-Or only knows Integer variables, setting variable to integer type";
      model_->setColumnIsInteger(index, true);
    }
    else
      model_->setColumnIsInteger(index, true);
#else
    glp_set_col_kind(lp_problem_, index + 1, (int)type);
#endif
  }

  LPWrapper::VariableType LPWrapper::getColumnType(Int index)
  {
#if COINOR_SOLVER == 1
    if (model_->isInteger(index))
    {
      return INTEGER;
    }
    else
      return CONTINUOUS;
#else
    return (VariableType)glp_get_col_kind(lp_problem_, index + 1);
#endif
  }

  void LPWrapper::setObjective(Int index, double obj_value)
  {
#if COINOR_SOLVER == 1
    model_->setObjective(index, obj_value);
#else
    glp_set_obj_coef(lp_problem_, index + 1, obj_value);
#endif
  }

  void LPWrapper::setObjectiveSense(LPWrapper::Sense sense)
  {
#if COINOR_SOLVER == 1
    if (sense == LPWrapper::MIN)
      model_->setOptimizationDirection(1);
    else
      model_->setOptimizationDirection(-1); // -1 maximize
#else
    glp_set_obj_dir(lp_problem_, (int)sense);
#endif
  }

  Int LPWrapper::getNumberOfColumns()
  {
#if COINOR_SOLVER == 1
    return model_->numberColumns();
#else
    return glp_get_num_cols(lp_problem_);
#endif
  }

  Int LPWrapper::getNumberOfRows()
  {
#if COINOR_SOLVER == 1
    return model_->numberRows();
#else
    return glp_get_num_rows(lp_problem_);
#endif
  }

  String LPWrapper::getColumnName(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->getColumnName(index);
#else
    return String(glp_get_col_name(lp_problem_, index + 1));
#endif
  }

  String LPWrapper::getRowName(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->getRowName(index);
#else
    return String(glp_get_row_name(lp_problem_, index + 1));
#endif
  }

  Int LPWrapper::getRowIndex(const String& name)
  {
#if COINOR_SOLVER == 1
    return model_->row(name.c_str());
#else
    glp_create_index(lp_problem_);
    return glp_find_row(lp_problem_, name.c_str()) - 1;
#endif
  }

  Int LPWrapper::getColumnIndex(const String& name)
  {
#if COINOR_SOLVER == 1
    return model_->column(name.c_str());
#else
    glp_create_index(lp_problem_);
    return glp_find_col(lp_problem_, name.c_str()) - 1;
#endif
  }

  LPWrapper::SOLVER LPWrapper::getSolver() const
  {
    return solver_;
  }

  void LPWrapper::readProblem(const String& filename, const String& format) // format=(LP,MPS,GLPK)
  {
#if COINOR_SOLVER == 1
    // delete old model and create a new model in its place (using same ptr)
    delete model_;
    model_ = new CoinModel(filename.c_str());
#else      
    // delete old model and create a new model in its place (using same ptr)
    glp_erase_prob(lp_problem_);

    if (format == "LP")
    {
      glp_read_lp(lp_problem_, nullptr, filename.c_str());
    }
    else if (format == "MPS")
    {
      glp_read_mps(lp_problem_, GLP_MPS_FILE, nullptr, filename.c_str());
    }
    else if (format == "GLPK")
    {
      glp_read_prob(lp_problem_, 0, filename.c_str());
    }
    else
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid LP format, allowed are LP, MPS, GLPK");

#endif
  }

  void LPWrapper::writeProblem(const String& filename, const WriteFormat format) const
  {
#if COINOR_SOLVER == 1
    if (format == FORMAT_MPS)
    {
      model_->writeMps(filename.c_str());
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid LP format, allowed is MPS");
    }
#else
    if (format == FORMAT_LP)
    {
      glp_write_lp(lp_problem_, nullptr, filename.c_str());
    }
    else if (format == FORMAT_MPS)
    {
      glp_write_mps(lp_problem_, GLP_MPS_FILE, nullptr, filename.c_str());
    }
    else if (format == FORMAT_GLPK)
    {
      glp_write_prob(lp_problem_, 0, filename.c_str());
    }
    else
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid LP format, allowed are LP, MPS, GLPK");
#endif
  }

  Int LPWrapper::solve(SolverParam& solver_param, const Size
#if COINOR_SOLVER == 1
                       verbose_level
#else
                       /* verbose_level */
#endif
                       )
  {
    OPENMS_LOG_INFO << "Using solver '" << (solver_ == LPWrapper::SOLVER_GLPK ? "glpk" : "coinor") << "' ...\n";
#if COINOR_SOLVER == 1
//Removed ifdef and OsiOslSolverInterface because Windows couldn't find it/both flags failed. For linux on the other hand the flags worked. But afaik we prefer CLP as solver anyway so no need to look for different solvers.
//#ifdef COIN_HAS_CLP
    OsiClpSolverInterface solver;
//#elif COIN_HAS_OSL
//      OsiOslSolverInterface solver;
//#endif
    solver.loadFromCoinModel(*model_);
    /* Now let MIP calculate a solution */
    // Pass to solver
    CbcModel model(solver);
    model.setObjSense(model_->optimizationDirection()); // -1 = maximize, 1=minimize
    model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

    // Output details
    model.messageHandler()->setLogLevel(verbose_level > 1 ? 2 : 0);
    model.solver()->messageHandler()->setLogLevel(verbose_level > 1 ? 1 : 0);

    //CglProbing generator1;
    //generator1.setUsingObjective(true);
    CglGomory generator2;
    generator2.setLimit(300);
    CglKnapsackCover generator3;
    CglOddHole generator4;
    generator4.setMinimumViolation(0.005);
    generator4.setMinimumViolationPer(0.00002);
    generator4.setMaximumEntries(200);
    CglClique generator5;
    generator5.setStarCliqueReport(false);
    generator5.setRowCliqueReport(false);
    //CglFlowCover flowGen;
    CglMixedIntegerRounding mixedGen;

    // Add in generators (you should prefer the ones used often and disable the others as they increase solution time)
    //model.addCutGenerator(&generator1,-1,"Probing");
    model.addCutGenerator(&generator2, -1, "Gomory");
    model.addCutGenerator(&generator3, -1, "Knapsack");
    //model.addCutGenerator(&generator4,-1,"OddHole"); // seg faults...
    model.addCutGenerator(&generator5, -10, "Clique");
    //model.addCutGenerator(&flowGen,-1,"FlowCover");
    model.addCutGenerator(&mixedGen, -1, "MixedIntegerRounding");

    // Heuristics
    CbcRounding heuristic1(model);
    model.addHeuristic(&heuristic1);
    CbcHeuristicLocal heuristic2(model);
    model.addHeuristic(&heuristic2);

    // set maximum allowed CPU time before forced stop (dangerous!)
    //model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*1);
      
    // Do initial solve to continuous
    model.initialSolve();


    // solve
    model.branchAndBound();
    // if (verbose_level > 0) OPENMS_LOG_INFO << " Branch and cut took " << CoinCpuTime()-time1 << " seconds, "
    //                                        << model.getNodeCount()<<" nodes with objective "
    //                                        << model.getObjValue()
    //                                        << (!model.status() ? " Finished" : " Not finished")
    //                                        << std::endl;
    for (Int i = 0; i < model_->numberColumns(); ++i)
    {
      solution_.push_back(model.solver()->getColSolution()[i]);
    }
    OPENMS_LOG_INFO << (model.isProvenOptimal() ? "Optimal solution found!" : "No solution found!") << "\n";
    return model.status();
#else

    glp_iocp solver_param_glp;
    glp_init_iocp(&solver_param_glp);

    solver_param_glp.msg_lev = solver_param.message_level;
    solver_param_glp.br_tech = solver_param.branching_tech;
    solver_param_glp.bt_tech = solver_param.backtrack_tech;
    solver_param_glp.pp_tech = solver_param.preprocessing_tech;
    if (solver_param.enable_feas_pump_heuristic)
      solver_param_glp.fp_heur = GLP_ON;
    if (solver_param.enable_gmi_cuts)
      solver_param_glp.gmi_cuts = GLP_ON;
    if (solver_param.enable_mir_cuts)
      solver_param_glp.mir_cuts = GLP_ON;
    if (solver_param.enable_cov_cuts)
      solver_param_glp.cov_cuts = GLP_ON;
    if (solver_param.enable_clq_cuts)
      solver_param_glp.clq_cuts = GLP_ON;
    solver_param_glp.mip_gap = solver_param.mip_gap;
    solver_param_glp.tm_lim = solver_param.time_limit;
    solver_param_glp.out_frq = solver_param.output_freq;
    solver_param_glp.out_dly = solver_param.output_delay;
    if (solver_param.enable_presolve)
      solver_param_glp.presolve = GLP_ON;
    if (solver_param.enable_binarization)
      solver_param_glp.binarize = GLP_ON; // only with presolve

    return glp_intopt(lp_problem_, &solver_param_glp);
#endif
  }

  LPWrapper::SolverStatus LPWrapper::getStatus()
  {
#if COINOR_SOLVER == 1
    return LPWrapper::UNDEFINED; // solver lokale Variable, braucht man diese Abfrage
#else
    Int status = glp_mip_status(lp_problem_);
    switch (status)
    {
    case 4:
      return LPWrapper::NO_FEASIBLE_SOL;

    case 5:
      return LPWrapper::OPTIMAL;

    case 2:
      return LPWrapper::FEASIBLE;

    default:
      return LPWrapper::UNDEFINED;
    }
#endif
  }

  double LPWrapper::getObjectiveValue()
  {
#if COINOR_SOLVER == 1
    double const * const obj = model_->objectiveArray();
    double obj_val = 0.;
    for (Int i = 0; i < model_->numberColumns(); ++i)
    {
      obj_val += obj[i] * getColumnValue(i);
    }
    return obj_val;
#else
    return glp_mip_obj_val(lp_problem_);
#endif
  }

  double LPWrapper::getColumnValue(Int index)
  {
#if COINOR_SOLVER == 1
    return solution_[index];
#else
    // glpk uses arrays beginning at pos 1, so we need to shift
    return glp_mip_col_val(lp_problem_, index + 1);
#endif
  }

  double LPWrapper::getColumnUpperBound(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->getColumnUpper(index);
#else
    return glp_get_col_ub(lp_problem_, index + 1);
#endif
  }

  double LPWrapper::getColumnLowerBound(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->getColumnLower(index);
#else
    return glp_get_col_lb(lp_problem_, index + 1);
#endif
  }

  double LPWrapper::getRowUpperBound(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->getRowUpper(index);
#else
    return glp_get_row_ub(lp_problem_, index + 1);
#endif
  }

  double LPWrapper::getRowLowerBound(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->getRowLower(index);
#else
    return glp_get_row_lb(lp_problem_, index + 1);
#endif
  }

  double LPWrapper::getObjective(Int index)
  {
#if COINOR_SOLVER == 1
    return model_->objective(index);
#else
    return glp_get_obj_coef(lp_problem_, index + 1);
#endif
  }

  LPWrapper::Sense LPWrapper::getObjectiveSense()
  {
#if COINOR_SOLVER == 1
    if (model_->optimizationDirection() == 1)
      return LPWrapper::MIN;
    else
      return LPWrapper::MAX;
#else
    if (glp_get_obj_dir(lp_problem_) == 1)
      return LPWrapper::MIN;
    else
      return LPWrapper::MAX;
#endif
  }

  Int LPWrapper::getNumberOfNonZeroEntriesInRow(Int idx)
  {
#if COINOR_SOLVER == 1
    const Int n_cols = getNumberOfColumns();
    std::vector<Int> ind(n_cols);
    std::vector<double> values(n_cols);
    Int nonzeroentries = 0;
    model_->getRow(idx, ind.data(), values.data());
    for (Int i = 0; i < n_cols; i++)
    {
      nonzeroentries += values[i] != 0 ? 1 : 0;
    }
    return nonzeroentries;
#else
    /* Non-zero coefficient count in the row. */
    // glpk uses arrays beginning at pos 1, so we need to shift
    return glp_get_mat_row(lp_problem_, idx + 1, nullptr, nullptr);
#endif
  }

  void LPWrapper::getMatrixRow(Int idx, std::vector<Int>& indexes)
  {
#if COINOR_SOLVER == 1
    indexes.clear();
    Int n_cols = getNumberOfColumns();
    std::vector<Int> ind(n_cols);
    std::vector<double> values(n_cols);
    model_->getRow(idx, ind.data(), values.data());
    for (Int i = 0; i < n_cols; ++i)
    {
      if (values[i] != 0)
        indexes.push_back(ind[i]);
    }
#else
    Int size = getNumberOfNonZeroEntriesInRow(idx);
    std::vector<int> ind(size + 1);
    glp_get_mat_row(lp_problem_, idx + 1, ind.data(), nullptr);
    indexes.clear();
    for (Int i = 1; i <= size; ++i)
    {
      indexes.push_back(ind[i] - 1);
    }
#endif
  }

} // namespace
