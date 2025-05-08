// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors:  Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ParameterInformation.h>

namespace OpenMS
{
  ParameterInformation::ParameterInformation(const String& n, ParameterTypes t, const String& arg, const ParamValue& def, const String& desc, bool req, bool adv, const StringList& tag_values) :
    name(n),
    type(t),
    default_value(def),
    description(desc),
    argument(arg),
    required(req),
    advanced(adv),
    tags(tag_values),
    valid_strings(),
    min_int(-std::numeric_limits<Int>::max()),
    max_int(std::numeric_limits<Int>::max()),
    min_float(-std::numeric_limits<double>::max()),
    max_float(std::numeric_limits<double>::max())
  {
  }

  ParameterInformation::ParameterInformation() :
    name(),
    type(NONE),
    default_value(),
    description(),
    argument(),
    required(true),
    advanced(false),
    tags(),
    valid_strings(),
    min_int(-std::numeric_limits<Int>::max()),
    max_int(std::numeric_limits<Int>::max()),
    min_float(-std::numeric_limits<double>::max()),
    max_float(std::numeric_limits<double>::max())
  {
  }

  ParameterInformation& ParameterInformation::operator=(const ParameterInformation& rhs)
  {
    if (&rhs == this)
    {
      return *this;
    }
    name = rhs.name;
    type = rhs.type;
    default_value = rhs.default_value;
    description = rhs.description;
    argument = rhs.argument;
    required = rhs.required;
    advanced = rhs.advanced;
    tags = rhs.tags;
    valid_strings = rhs.valid_strings;
    min_int = rhs.min_int;
    max_int = rhs.max_int;
    min_float = rhs.min_float;
    max_float = rhs.max_float;

    return *this;
  }

}
