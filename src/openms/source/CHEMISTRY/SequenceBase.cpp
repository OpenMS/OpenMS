// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SequenceBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS
{

  template<>
  double SequenceBase<AASequence>::getMonoWeight(Int charge) const
  {
    return derived().getMonoWeight(Residue::ResidueType::Full, charge);
  }

  template<>
  double SequenceBase<NASequence>::getMonoWeight(Int charge) const
  {
    return derived().getMonoWeight(NASequence::NASFragmentType::Full, charge);
  }

  template<>
  double SequenceBase<AASequence>::getAverageWeight(Int charge) const
  {
    return derived().getAverageWeight(Residue::ResidueType::Full, charge);
  }

  template<>
  double SequenceBase<NASequence>::getAverageWeight(Int charge) const
  {
    return derived().getAverageWeight(NASequence::NASFragmentType::Full, charge);
  }

  template<>
  EmpiricalFormula SequenceBase<AASequence>::getFormula(Int charge) const
  {
    return derived().getFormula(Residue::ResidueType::Full, charge);
  }

  template<>
  EmpiricalFormula SequenceBase<NASequence>::getFormula(Int charge) const
  {
    return derived().getFormula(NASequence::NASFragmentType::Full, charge);
  }

  // Explicit template instantiations after specializations
  template class SequenceBase<AASequence>;
  template class SequenceBase<NASequence>;

} // namespace OpenMS