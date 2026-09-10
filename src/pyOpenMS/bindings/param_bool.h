// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Müller $
// $Authors: Tom David Müller $
// --------------------------------------------------------------------------

/**
 * @file param_bool.h
 * @brief Python bool <-> OpenMS boolean parameter translation for the Param bindings
 *
 * OpenMS has no boolean ParamValue. A boolean parameter ("flag") is a
 * STRING_VALUE holding "true"/"false"; the C++ side reads it by content
 * (ParamValue::toBool(), `getValue(...) == "true"`) and only uses the entry's
 * restrictions (valid_strings {"true","false"}) when rendering INI/CTD/CWL files.
 *
 * pyOpenMS follows the same rule so that a boolean parameter is a Python bool
 * on every path (getValue, [], items, ParamEntry.value, ...), no matter how the
 * Param was built (algorithm defaults, INI file, setDefaults/merge/insert,
 * user assignment of True or of the string "true"):
 *
 *   read : STRING_VALUE exactly "true"/"false" -> bool, unless the entry's
 *          restrictions contain other values (tri-state parameters such as
 *          "auto,true,false" stay str on every path)
 *   write: bool -> "true"/"false"; an entry without restrictions additionally
 *          gets valid_strings {"true","false"} so that the C++ file writers
 *          recognise it as a flag
 *   restrictions: {"true","false"} in either order are shown as [True, False]
 *   type : Param.getValueType() reports ValueType.BOOL_VALUE (Python-only)
 */

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>

#include <algorithm>
#include <string>
#include <vector>

#include "type_casters/openms_datavalue_caster.h"

namespace pyopenms_param_bool {

namespace nb = nanobind;

/// Signature string for values accepted by the Param write paths (kept in sync
/// with the ParamValue type caster in type_casters/openms_datavalue_caster.h).
constexpr const char* kParamValueSig = "None | bool | int | float | str | bytes | list[str] | list[int] | list[float]";

/// Python-only ValueType member reported for boolean parameters. Outside the
/// C++ range (0-6); ValueType has the fixed underlying type unsigned char and
/// never enters C++ from Python, so the value cannot reach OpenMS.
constexpr OpenMS::ParamValue::ValueType kParamBoolValueType = static_cast<OpenMS::ParamValue::ValueType>(7);

inline const std::vector<std::string>& paramBoolValidStrings()
{
    static const std::vector<std::string> valid{"true", "false"};
    return valid;
}

inline bool isBoolString(const std::string& s)
{
    return s == "true" || s == "false";
}

/// True if the restriction list is exactly {"true","false"} in either order.
inline bool paramValidStringsAreBool(const std::vector<std::string>& valid_strings)
{
    return valid_strings.size() == 2
        && isBoolString(valid_strings[0]) && isBoolString(valid_strings[1])
        && valid_strings[0] != valid_strings[1];
}

/// True if the entry may hold a boolean: a string parameter whose restrictions
/// (if any) allow nothing but "true"/"false".
inline bool paramEntryIsBoolTyped(const OpenMS::Param::ParamEntry& entry)
{
    return entry.value.valueType() == OpenMS::ParamValue::STRING_VALUE
        && std::all_of(entry.valid_strings.begin(), entry.valid_strings.end(), isBoolString);
}

/// The single predicate behind bool reads and Param.getValueType(): the entry
/// is bool-typed and currently holds "true" or "false". Any other string (e.g.
/// an invalid "maybe" that checkDefaults is meant to reject) stays a str.
inline bool paramEntryReadsAsBool(const OpenMS::Param::ParamEntry& entry)
{
    return paramEntryIsBoolTyped(entry) && isBoolString(static_cast<std::string>(entry.value));
}

/// Read side: Python True/False for a boolean parameter, otherwise the plain caster result.
inline nb::object paramEntryValueToPython(const OpenMS::Param::ParamEntry& entry)
{
    if (paramEntryReadsAsBool(entry))
    {
        return nb::bool_(static_cast<std::string>(entry.value) == "true");
    }
    return nb::cast(entry.value);
}

/// Restrictions as seen from Python: [True, False] for a boolean restriction, else list[str].
inline nb::object paramValidStringsToPython(const std::vector<std::string>& valid_strings)
{
    if (paramValidStringsAreBool(valid_strings))
    {
        nb::list result;
        result.append(nb::bool_(true));
        result.append(nb::bool_(false));
        return result;
    }
    return nb::cast(valid_strings);
}

/// Python sequence -> valid strings. bool elements become "true"/"false" and a
/// list consisting of exactly {True, False} is stored in the canonical order
/// {"true","false"} that the C++ file writers recognise. str/bytes are kept.
inline std::vector<std::string> paramValidStringsFromPython(nb::handle strings)
{
    if (PyUnicode_Check(strings.ptr()) || PyBytes_Check(strings.ptr()))
    {
        throw nb::type_error("Valid strings must be a list of str/bytes/bool, not a single string");
    }
    std::vector<std::string> result;
    for (nb::handle item : strings)
    {
        if (PyBool_Check(item.ptr()))
        {
            result.emplace_back(item.ptr() == Py_True ? "true" : "false");
            continue;
        }
        std::string s;
        if (!nb::try_cast(item, s))
        {
            const std::string type_name = nb::cast<std::string>(nb::str(item.type().attr("__name__")));
            throw nb::type_error(("Valid strings must be str, bytes or bool, got " + type_name).c_str());
        }
        result.push_back(std::move(s));
    }
    if (paramValidStringsAreBool(result))
    {
        result = paramBoolValidStrings();
    }
    return result;
}

/// Write side: convert a Python object to ParamValue, remembering whether it
/// was a bool (nanobind overload dispatch on `const ParamValue&` would lose
/// that information, and the caster cannot touch the entry's restrictions).
inline OpenMS::ParamValue paramValueFromPython(nb::handle value, bool& was_bool)
{
    was_bool = PyBool_Check(value.ptr());
    OpenMS::ParamValue pv;
    if (!nb::try_cast(value, pv))
    {
        const std::string type_name = nb::cast<std::string>(nb::str(value.type().attr("__name__")));
        throw nb::type_error(("Param value must be " + std::string(kParamValueSig) + ", got " + type_name).c_str());
    }
    return pv;
}

/// Single entry point for every Param write path. Sets the value and, when the
/// source was a Python bool and the entry carries no restriction yet, adds the
/// {"true","false"} restriction so that INI/CTD writers render it as a flag.
/// Existing restrictions are never overwritten (Param::setValue keeps the
/// valid_strings of an existing entry).
inline void paramSetValueFromPython(OpenMS::Param& param, const std::string& key, nb::handle value,
                                    const std::string& description = "",
                                    const std::vector<std::string>& tags = std::vector<std::string>())
{
    bool was_bool = false;
    const OpenMS::ParamValue pv = paramValueFromPython(value, was_bool);
    param.setValue(key, pv, description, tags);
    if (was_bool && param.getEntry(key).valid_strings.empty())
    {
        // cannot throw: the entry is a STRING_VALUE and the strings contain no comma
        param.setValidStrings(key, paramBoolValidStrings());
    }
}

} // namespace pyopenms_param_bool
