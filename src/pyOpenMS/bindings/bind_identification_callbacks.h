// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

// What the bindings of identification data and of the functions that call Python back with its records share
// (the identification data module binds RunView; e.g. the IDFilter predicates get runs as RunViews)

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <nanobind/nanobind.h>
#include <string>

namespace pyopenms_identification
{
namespace nb = nanobind;
using ID = OpenMS::IdentificationData;

/// Python truthiness of a callback result (accepts numpy.bool_, 0/1, None, ...).
inline bool truthy(const nb::object& value)
{
  const int result = PyObject_IsTrue(value.ptr());
  if (result < 0) throw nb::python_error();
  return result != 0;
}

/// A run inside a dataset, addressed by UUID: it cannot dangle when the dataset's runs change,
/// and every edit goes to the live run, so IDs are only ever allocated in one place.
struct RunView
{
  nb::object owner; // the IdentificationData (an owned object or a view into a feature/consensus map)
  std::string uuid;
  ID::Run& run() const
  {
    auto* found = nb::cast<ID&>(owner).findRunByUuid(uuid);
    if (! found) throw nb::key_error(("The run " + uuid + " is no longer part of the dataset").c_str());
    return *found;
  }
};
} // namespace pyopenms_identification
