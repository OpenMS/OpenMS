// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/config.h>

#include <Highs.h>

namespace OpenMS
{
  LPWrapper::LPWrapper()
  {
    highs_ = std::make_unique<Highs>();
    highs_->setOptionValue("output_flag", false);
    // Add an empty LP model
    HighsLp lp;
    lp.num_col_ = 0;
    lp.num_row_ = 0;
    highs_->passModel(std::move(lp));
  }

  LPWrapper::~LPWrapper()
  {
    // unique_ptr handles cleanup
  }

  Int LPWrapper::addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values, const std::string& name) // return index
  {
    if (row_indices.size() != row_values.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Indices and values vectors differ in size");
    }
    // Add a row with no bounds (free row) by default
    std::vector<HighsInt> indices(row_indices.begin(), row_indices.end());
    highs_->addRow(-kHighsInf, kHighsInf, (HighsInt)indices.size(), indices.data(), row_values.data());
    const Int index = highs_->getNumRow() - 1;
    highs_->passRowName(index, name);
    return index;
  }

  Int LPWrapper::addColumn()
  {
    // Add a column with bounds fixed at 0
    highs_->addCol(0.0, 0.0, 0.0, 0, nullptr, nullptr);
    return highs_->getNumCol() - 1;
  }

  Int LPWrapper::addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values, const std::string& name)
  {
    if (column_indices.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Column indices for Row are empty");
    }
    if (column_indices.size() != column_values.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Indices and values vectors differ in size");
    }
    std::vector<HighsInt> indices(column_indices.begin(), column_indices.end());
    highs_->addCol(0.0, -kHighsInf, kHighsInf, (HighsInt)indices.size(), indices.data(), column_values.data());
    const Int index = highs_->getNumCol() - 1;
    highs_->passColName(index, name);
    return index;
  }

  Int LPWrapper::addRow(const std::vector<Int>& row_indices, const std::vector<double>& row_values,
                        const std::string& name, double lower_bound, double upper_bound, Type type)
  {
    const Int index = addRow(row_indices, row_values, name);
    setRowBounds(index, lower_bound, upper_bound, type);
    return index; // in addRow index is decreased already
  }

  Int LPWrapper::addColumn(const std::vector<Int>& column_indices, const std::vector<double>& column_values,
                           const std::string& name, double lower_bound, double upper_bound, Type type) //return index
  {
    const Int index = addColumn(column_indices, column_values, name);
    setColumnBounds(index, lower_bound, upper_bound, type);
    return index; // in addColumn index is decreased already
  }

  void LPWrapper::deleteRow(Int index)
  {
    HighsInt idx = index;
    highs_->deleteRows(1, &idx);
  }

  void LPWrapper::setElement(Int row_index, Int column_index, double value)
  {
    if (row_index >= getNumberOfRows() || column_index >= getNumberOfColumns())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid index given",std::string("invalid column_index or row_index"));
    }
    highs_->changeCoeff(row_index, column_index, value);
  }

  double LPWrapper::getElement(Int row_index, Int column_index)
  {
    if (row_index >= getNumberOfRows() || column_index >= getNumberOfColumns())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid index given",std::string("invalid column_index or row_index"));
    }
  double value = 0.0;
  highs_->getCoeff(row_index, column_index, value);
  return value;
  }

  void LPWrapper::setColumnName(Int index, const std::string& name)
  {
    highs_->passColName(index, name);
  }

  void LPWrapper::setRowName(Int index, const std::string& name)
  {
    highs_->passRowName(index, name);
  }

  void LPWrapper::setColumnBounds(Int index, double lower_bound, double upper_bound, Type type)
  {
    switch (type)
    {
    case UNBOUNDED:
      highs_->changeColBounds(index, -kHighsInf, kHighsInf);
      break;
    case LOWER_BOUND_ONLY:
      highs_->changeColBounds(index, lower_bound, kHighsInf);
      break;
    case UPPER_BOUND_ONLY:
      highs_->changeColBounds(index, -kHighsInf, upper_bound);
      break;
    default: // double-bounded or fixed
      highs_->changeColBounds(index, lower_bound, upper_bound);
      break;
    }
  }

  void LPWrapper::setRowBounds(Int index, double lower_bound, double upper_bound, Type type)
  {
    switch (type)
    {
    case UNBOUNDED:
      highs_->changeRowBounds(index, -kHighsInf, kHighsInf);
      break;
    case LOWER_BOUND_ONLY:
      highs_->changeRowBounds(index, lower_bound, kHighsInf);
      break;
    case UPPER_BOUND_ONLY:
      highs_->changeRowBounds(index, -kHighsInf, upper_bound);
      break;
    default: // double-bounded or fixed
      highs_->changeRowBounds(index, lower_bound, upper_bound);
      break;
    }
  }

  void LPWrapper::setColumnType(Int index, VariableType type) // 1- continuous, 2- integer, 3- binary
  {
    if (index < 0 || index >= highs_->getNumCol())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, highs_->getNumCol());
    }
    // stage the change; applied in one batch by flushIntegrality_()
    if (staged_integrality_.empty())
    {
      const HighsLp& lp = highs_->getLp();
      for (HighsVarType t : lp.integrality_) staged_integrality_.push_back(static_cast<uint8_t>(t));
    }
    staged_integrality_.resize(highs_->getNumCol(), static_cast<uint8_t>(HighsVarType::kContinuous));
    staged_integrality_[index] = static_cast<uint8_t>(type == CONTINUOUS ? HighsVarType::kContinuous : HighsVarType::kInteger);
    if (type == BINARY)
    {
      highs_->changeColBounds(index, 0.0, 1.0);
    }
  }

  LPWrapper::VariableType LPWrapper::getColumnType(Int index)
  {
    // Get integrality from the model's integrality vector
    flushIntegrality_();
    const HighsLp& lp = highs_->getLp();
    if (lp.integrality_.empty() || lp.integrality_[index] == HighsVarType::kContinuous)
      return CONTINUOUS;
    else
      return INTEGER;
  }

  void LPWrapper::setObjective(Int index, double obj_value)
  {
  highs_->changeColCost(index, obj_value);
  }

  void LPWrapper::setObjectiveSense(LPWrapper::Sense sense)
  {
    if (sense == LPWrapper::MIN)
      highs_->changeObjectiveSense(ObjSense::kMinimize);
    else
      highs_->changeObjectiveSense(ObjSense::kMaximize);
  }

  Int LPWrapper::getNumberOfColumns()
  {
    return highs_->getNumCol();
  }

  Int LPWrapper::getNumberOfRows()
  {
    return highs_->getNumRow();
  }

  std::string LPWrapper::getColumnName(Int index)
  {
    std::string name;
    highs_->getColName(index, name);
    return StringUtils::toStr(name);
  }

  std::string LPWrapper::getRowName(Int index)
  {
    std::string name;
    highs_->getRowName(index, name);
    return StringUtils::toStr(name);
  }

  Int LPWrapper::getRowIndex(const std::string& name)
  {
    // Linear search through row names
    const HighsInt num_rows = highs_->getNumRow();
    for (HighsInt i = 0; i < num_rows; ++i)
    {
      std::string rname;
      highs_->getRowName(i, rname);
      if (rname == std::string(name.c_str()))
        return i;
    }
    return -1;
  }

  Int LPWrapper::getColumnIndex(const std::string& name)
  {
    // Linear search through column names
    const HighsInt num_cols = highs_->getNumCol();
    for (HighsInt i = 0; i < num_cols; ++i)
    {
      std::string cname;
      highs_->getColName(i, cname);
      if (cname == std::string(name.c_str()))
        return i;
    }
    return -1;
  }

  LPWrapper::SOLVER LPWrapper::getSolver() const
  {
    return SOLVER_HIGHS;
  }

  void LPWrapper::readProblem(const std::string& filename, const std::string& /*format*/)
  {
    // Clear and read a new model
    staged_integrality_.clear();
    highs_->clear();
    highs_->setOptionValue("output_flag", false);
    highs_->readModel(filename.c_str());
  }

  void LPWrapper::writeProblem(const std::string& filename, const WriteFormat format) const
  {
    (void)format;
    flushIntegrality_();
    // HiGHS can write MPS and LP formats; it auto-detects from extension
    highs_->writeModel(filename.c_str());
  }

  Int LPWrapper::solve(SolverParam& solver_param, const Size verbose_level)
  {
    OPENMS_LOG_INFO << "Using solver 'highs' ...\n";
    // Configure HiGHS options based on solver_param
    flushIntegrality_();
    highs_->setOptionValue("output_flag", verbose_level > 1);
    // set (or reset) on every call, so a limit from a previous solve() does not carry over; the limit applies per
    // solve(), since HiGHS' run clock accumulates over run() calls and is zeroed below
    highs_->setOptionValue("time_limit", solver_param.time_limit < (std::numeric_limits<Int>::max)()
                                           ? solver_param.time_limit / 1000.0 // convert ms to seconds
                                           : kHighsInf);
    // always pass the gap: HiGHS' own default (1e-4) would otherwise stop before proven optimality for mip_gap = 0
    highs_->setOptionValue("mip_rel_gap", solver_param.mip_gap);
    highs_->setOptionValue("presolve", solver_param.enable_presolve ? std::string("on") : std::string("off"));
    // Skip the probing and enumeration presolve rules. On the large binary set-packing ILPs OpenMS builds
    // (e.g. decharging), they spend most of the runtime building clique tables (50k-node instance:
    // 70 s -> 5 s, same optimum). Only switching off both helps; the remaining presolve rules stay on.
    HighsInt presolve_rule_off = 1 << kPresolveRuleProbing;
#if HIGHS_VERSION_MAJOR > 1 || HIGHS_VERSION_MINOR >= 13 // the enumeration rule was added in HiGHS 1.13
    presolve_rule_off |= 1 << kPresolveRuleEnumeration;
#endif
    highs_->setOptionValue("presolve_rule_off", presolve_rule_off);

    // Solve
    highs_->zeroAllClocks();
    HighsStatus status = highs_->run();

    // Store solution
    solution_.clear();
    if (status == HighsStatus::kOk || status == HighsStatus::kWarning)
    {
      const HighsSolution& sol = highs_->getSolution();
      solution_ = sol.col_value;
      HighsModelStatus model_status = highs_->getModelStatus();
      OPENMS_LOG_INFO << (model_status == HighsModelStatus::kOptimal ? "Optimal solution found!" : "Solution status: non-optimal") << "\n";
    }
    else
    {
      OPENMS_LOG_INFO << "No solution found!\n";
    }
    const Int raw_solver_return = static_cast<Int>(status);
    // make solver failures visible: without a feasible solution all column values read as zero,
    // which is indistinguishable from a genuine all-zero optimum (see issue #9944)
    const SolverStatus solution_status = getStatus();
    if (solution_status != LPWrapper::OPTIMAL && solution_status != LPWrapper::FEASIBLE)
    {
      OPENMS_LOG_WARN << "LPWrapper::solve(): the solver did not find a feasible solution (raw solver return code: " << raw_solver_return
                      << "). Solution values are not meaningful; check getStatus() before using them." << std::endl;
    }
    return raw_solver_return;
  }

  LPWrapper::SolverStatus LPWrapper::getStatus()
  {
    HighsModelStatus model_status = highs_->getModelStatus();
    switch (model_status)
    {
    case HighsModelStatus::kOptimal:
      return LPWrapper::OPTIMAL;
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
      return LPWrapper::FEASIBLE;
    case HighsModelStatus::kInfeasible:
    case HighsModelStatus::kUnbounded:
      return LPWrapper::NO_FEASIBLE_SOL;
    default:
      // e.g. a time/iteration limit was hit: an integer feasible (incumbent) solution may still exist
      if (highs_->getHighsInfo().primal_solution_status == kSolutionStatusFeasible)
      {
        return LPWrapper::FEASIBLE;
      }
      return LPWrapper::UNDEFINED;
    }
  }

  double LPWrapper::getObjectiveValue()
  {
    const HighsInfo& info = highs_->getHighsInfo();
    return info.objective_function_value;
  }

  double LPWrapper::getColumnValue(Int index)
  {
    if (index >= 0 && index < (Int)solution_.size())
      return solution_[index];
    return 0.0;
  }

  double LPWrapper::getColumnUpperBound(Int index)
  {
    const HighsLp& lp = highs_->getLp();
    return lp.col_upper_[index];
  }

  double LPWrapper::getColumnLowerBound(Int index)
  {
    const HighsLp& lp = highs_->getLp();
    return lp.col_lower_[index];
  }

  double LPWrapper::getRowUpperBound(Int index)
  {
    const HighsLp& lp = highs_->getLp();
    return lp.row_upper_[index];
  }

  double LPWrapper::getRowLowerBound(Int index)
  {
    const HighsLp& lp = highs_->getLp();
    return lp.row_lower_[index];
  }

  double LPWrapper::getObjective(Int index)
  {
    const HighsLp& lp = highs_->getLp();
    return lp.col_cost_[index];
  }

  LPWrapper::Sense LPWrapper::getObjectiveSense()
  {
    const HighsLp& lp = highs_->getLp();
    if (lp.sense_ == ObjSense::kMinimize)
      return LPWrapper::MIN;
    else
      return LPWrapper::MAX;
  }

  Int LPWrapper::getNumberOfNonZeroEntriesInRow(Int idx)
  {
    HighsInt num_row;
    HighsInt num_nz;
  HighsInt row_start;
  highs_->getRows(idx, idx, num_row, nullptr, nullptr, num_nz, &row_start, nullptr, nullptr);
    return num_nz;
  }

  void LPWrapper::getMatrixRow(Int idx, std::vector<Int>& indexes)
  {
    indexes.clear();
    HighsInt num_row;
    HighsInt num_nz;
  HighsInt row_start;
  highs_->getRows(idx, idx, num_row, nullptr, nullptr, num_nz, &row_start, nullptr, nullptr);
    if (num_nz == 0) return;
    std::vector<double> row_lower(1), row_upper(1);
  std::vector<HighsInt> starts(1);
    std::vector<HighsInt> ind(num_nz);
    std::vector<double> values(num_nz);
    HighsInt num_nz_out;
  highs_->getRows(idx, idx, num_row, row_lower.data(), row_upper.data(), num_nz_out, starts.data(), ind.data(), values.data());
    for (HighsInt i = 0; i < num_nz_out; ++i)
    {
      indexes.push_back(static_cast<Int>(ind[i]));
    }
  }

  void LPWrapper::flushIntegrality_() const
  {
    if (staged_integrality_.empty())
    {
      return;
    }
    const HighsInt num_col = highs_->getNumCol();
    // columns added after the last setColumnType() call are continuous
    staged_integrality_.resize(num_col, static_cast<uint8_t>(HighsVarType::kContinuous));
    std::vector<HighsVarType> integrality(num_col);
    for (HighsInt i = 0; i < num_col; ++i)
    {
      integrality[i] = static_cast<HighsVarType>(staged_integrality_[i]);
    }
    if (num_col > 0)
    {
      highs_->changeColsIntegrality(0, num_col - 1, integrality.data());
    }
    staged_integrality_.clear();
  }

} // namespace
