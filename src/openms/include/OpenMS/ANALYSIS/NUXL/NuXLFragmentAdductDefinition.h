// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <vector>
#include <map>
#include <set>

namespace OpenMS
{  
struct OPENMS_DLLAPI NuXLFragmentAdductDefinition
{
  EmpiricalFormula formula; // formula
  String name;  // name used in annotation
  double mass = 0;

  NuXLFragmentAdductDefinition() = default;

  NuXLFragmentAdductDefinition(const NuXLFragmentAdductDefinition&) = default;

  NuXLFragmentAdductDefinition(NuXLFragmentAdductDefinition&&) = default;

  NuXLFragmentAdductDefinition& operator=(const NuXLFragmentAdductDefinition&) = default;

  NuXLFragmentAdductDefinition& operator=(NuXLFragmentAdductDefinition&&) = default;

  bool operator<(const NuXLFragmentAdductDefinition& other) const;

  bool operator==(const NuXLFragmentAdductDefinition& other) const;
};

}

