// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "OpenMS/KERNEL/Feature.h"
#include "OpenMS/METADATA/PeptideHit.h"

#include <optional>

namespace OpenMS::PipEcho::Util
{

/******************************************************************************/
/**
 * Return the first peptide hit from a feature.
 *
 * This code is here to isolate checking hits because the process
 * might change in a future version of OpenMS.
 */
std::optional<PeptideHit> feature_hit(const Feature&);

/******************************************************************************/
/**
 * Compute the mass error in PPM for the given feature.
 */
std::optional<double> feature_mass_error(const Feature&);

/******************************************************************************/
/**
 * Function object that wraps around another function object `F`.
 * When invoked, this function object will dereference a pointer and
 * invoke `F` with the result.
 *
 */
template<typename F>
struct dref_fn_t
{
  F f;

  /// operator() :: (T -> R) -> T* -> R
  template<typename T>
  constexpr decltype(auto) operator()(T&& t)
  {
    return std::invoke(f, *t);
  };
};

/******************************************************************************/
/**
 * See `dref_fn_t`.
 */
template<typename F>
constexpr dref_fn_t<std::decay_t<F>> dref_fn(F&& f)
{
  return {std::forward<F>(f)};
}

} // namespace OpenMS::PipEcho::Util
