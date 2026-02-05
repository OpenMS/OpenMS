/**
 * @file all_casters.h
 * @brief Master include file for all pyOpenMS type casters
 *
 * This file includes all custom type casters required for converting
 * between Python and OpenMS C++ types.
 *
 * Include this file before any nanobind class definitions.
 */

#pragma once

#include <nanobind/nanobind.h>
// Use custom std::string caster that accepts both str and bytes
#include "std_string_bytes_caster.h"
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/ndarray.h>

// OpenMS-specific type casters (prefixed to avoid collision with system headers)
#include "openms_string_caster.h"
#include "openms_dposition_caster.h"
#include "openms_datavalue_caster.h"
#include "openms_stl_caster.h"

namespace nb = nanobind;
