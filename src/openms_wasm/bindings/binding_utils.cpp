// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "binding_utils.h"

#include <cmath>
#include <limits>

namespace OpenMS::WASM
{
  namespace
  {
    /// Whether the number @p d is integral and fits into an int (the integers of a Param)
    bool fitsInt(double d)
    {
      return d == std::floor(d) && d >= std::numeric_limits<int>::min() && d <= std::numeric_limits<int>::max();
    }

    /// The element type of a non-empty JavaScript array: a list of strings, ints or doubles.
    ParamValue::ValueType listTypeOf(const val& array)
    {
      bool all_integral = true;
      const unsigned n = array["length"].as<unsigned>();
      for (unsigned i = 0; i < n; ++i)
      {
        const val element = array[i];
        if (element.isString())
        {
          return ParamValue::STRING_LIST;
        }
        all_integral = all_integral && element.isNumber() && fitsInt(element.as<double>());
      }
      return all_integral ? ParamValue::INT_LIST : ParamValue::DOUBLE_LIST;
    }
  } // namespace

  val toJS(const ParamValue& v, bool as_bool)
  {
    switch (v.valueType())
    {
      case ParamValue::STRING_VALUE:
      {
        const std::string s = v.toString();
        if (as_bool && (s == "true" || s == "false"))
        {
          return val(s == "true");
        }
        return val(s);
      }
      case ParamValue::INT_VALUE:
        return val(static_cast<double>(static_cast<int>(v)));
      case ParamValue::DOUBLE_VALUE:
        return val(static_cast<double>(v));
      case ParamValue::STRING_LIST:
        return toArray(v.toStringVector());
      case ParamValue::INT_LIST:
        return toArray(v.toIntVector());
      case ParamValue::DOUBLE_LIST:
        return toArray(v.toDoubleVector());
      default:
        return val::null();
    }
  }

  val toJS(const DataValue& v)
  {
    switch (v.valueType())
    {
      case DataValue::STRING_VALUE:
        return val(v.toString());
      case DataValue::INT_VALUE:
        // an integer meta value may exceed 2^31 (e.g. 64-bit ids): keep it as double
        return val(static_cast<double>(static_cast<long long>(v)));
      case DataValue::DOUBLE_VALUE:
        return val(static_cast<double>(v));
      case DataValue::STRING_LIST:
        return toArray(v.toStringList());
      case DataValue::INT_LIST:
        return toArray(v.toIntList());
      case DataValue::DOUBLE_LIST:
        return toArray(v.toDoubleList());
      default:
        return val::null();
    }
  }

  ParamValue paramValueFromJS(const val& v, ParamValue::ValueType type)
  {
    if (isNullish(v))
    {
      return ParamValue();
    }
    if (type == ParamValue::EMPTY_VALUE)
    {
      if (v.isString() || v.isTrue() || v.isFalse())
      {
        type = ParamValue::STRING_VALUE;
      }
      else if (v.isNumber())
      {
        type = fitsInt(v.as<double>()) ? ParamValue::INT_VALUE : ParamValue::DOUBLE_VALUE;
      }
      else if (val::global("Array").call<bool>("isArray", v))
      {
        type = listTypeOf(v);
      }
      else
      {
        throwJSError("Unsupported parameter value; expected a number, string, boolean or array.");
      }
    }

    switch (type)
    {
      case ParamValue::STRING_VALUE:
        if (v.isTrue() || v.isFalse())
        {
          return ParamValue(v.isTrue() ? "true" : "false");
        }
        if (v.isNumber())
        {
          return ParamValue(val::global("String")(v).as<std::string>());
        }
        return ParamValue(v.as<std::string>());
      case ParamValue::INT_VALUE:
      {
        const double d = v.as<double>();
        if (!fitsInt(d))
        {
          throwJSError("Integer parameter value expected (32 bit), got " + val::global("String")(v).as<std::string>());
        }
        return ParamValue(static_cast<int>(d));
      }
      case ParamValue::DOUBLE_VALUE:
        return ParamValue(v.as<double>());
      case ParamValue::STRING_LIST:
        return ParamValue(stringsFromJS(v));
      case ParamValue::INT_LIST:
        return ParamValue(numbersFromJS<int>(v));
      case ParamValue::DOUBLE_LIST:
        return ParamValue(numbersFromJS<double>(v));
      default:
        return ParamValue();
    }
  }

  DataValue dataValueFromJS(const val& v)
  {
    // A meta value holds 64-bit integers, a parameter only 32-bit ones: take integral numbers
    // up to 2^53 (exact in JavaScript) directly
    if (v.isNumber())
    {
      const double d = v.as<double>();
      if (d == std::floor(d) && std::fabs(d) <= 9007199254740992.0)
      {
        return DataValue(static_cast<long long>(d));
      }
    }
    return DataValue(paramValueFromJS(v));
  }

  JS::ParamValues paramToObject(const Param& p)
  {
    val result = val::object();
    for (auto it = p.begin(); it != p.end(); ++it)
    {
      result.set(it.getName(), toJS(it->value, it->isBool()));
    }
    return result;
  }

  void updateParamFromObject(Param& p, const JS::ParamValues& values)
  {
    const val keys = val::global("Object").call<val>("keys", values);
    const unsigned n = keys["length"].as<unsigned>();
    for (unsigned i = 0; i < n; ++i)
    {
      const std::string key = keys[i].as<std::string>();
      if (!p.exists(key))
      {
        throwJSError("Unknown parameter '" + key + "'");
      }
      const Param::ParamEntry& entry = p.getEntry(key);
      p.setValue(key, paramValueFromJS(values[key], entry.value.valueType()), entry.description,
                 std::vector<std::string>(entry.tags.begin(), entry.tags.end()));
    }
  }
} // namespace OpenMS::WASM
