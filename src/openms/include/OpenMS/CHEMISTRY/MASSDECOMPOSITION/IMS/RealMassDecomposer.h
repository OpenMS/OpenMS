// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <utility>
#include <map>
#include <memory>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>

namespace OpenMS
{
  namespace ims
  {

    /**
      @brief Decomposes a real-valued mass over a real-valued alphabet
             of @ref Weights, within a configurable absolute mass error.

      Translates the real-valued problem into integer mass decomposition
      against the @ref Weights' scaled integer alphabet: it computes the
      integer-mass window that brackets the user-supplied
      @c [mass - error, mass + error] range (accounting for the
      precision-induced rounding errors stored on @ref Weights), runs an
      @ref IntegerMassDecomposer over each integer mass in the window,
      then drops every integer decomposition whose reconstructed real
      parent mass still differs from @p mass by more than @p error.

      @ingroup Analysis_DeNovo
    */
    class OPENMS_DLLAPI RealMassDecomposer
    {
public:

      /// Underlying integer decomposer.
      typedef IntegerMassDecomposer<> integer_decomposer_type;

      /// Integer mass values consumed by the integer decomposer.
      typedef integer_decomposer_type::value_type integer_value_type;

      /// Result type of the integer decomposer (collection of decompositions).
      typedef integer_decomposer_type::decompositions_type decompositions_type;

      /// Counter type returned by @ref getNumberOfDecompositions.
      typedef unsigned long long number_of_decompositions_type;

      /// Per-alphabet-entry count constraint: alphabet index -> (min count, max count); both bounds inclusive.
      typedef std::map<unsigned int, std::pair<unsigned int, unsigned int> > constraints_type;

      /**
        @brief Construct from an alphabet of @ref Weights.

        The alphabet's precision and rounding-error bounds are captured
        at construction time and reused on every subsequent
        @ref getDecompositions / @ref getNumberOfDecompositions call.

        @param[in] weights Real-valued alphabet to decompose against.
      */
      explicit RealMassDecomposer(const Weights & weights);

      /**
        @brief All decompositions of @p mass within absolute tolerance @p error.

        Returns every alphabet-multiplicity vector whose reconstructed
        parent mass lies in @c [mass - error, mass + error].

        @param[in] mass  Target real-valued mass.
        @param[in] error Absolute mass tolerance.
        @return Decompositions whose reconstructed real parent mass is
                within @p error of @p mass.
      */
      decompositions_type getDecompositions(double mass, double error);

      /**
        @brief Same as @ref getDecompositions(double,double), additionally filtered by per-alphabet count bounds.

        For each entry @c (idx, (lo, hi)) in @p constraints, only
        decompositions @c d with @c lo @c <= @c d[idx] @c <= @c hi are
        kept. Constraints whose alphabet index does not appear in @p
        decompositions are simply skipped.

        @param[in] mass        Target real-valued mass.
        @param[in] error       Absolute mass tolerance.
        @param[in] constraints Per-alphabet-entry count bounds.
        @return Decompositions that satisfy both the mass tolerance and
                every supplied count constraint.
      */
      decompositions_type getDecompositions(double mass, double error, const constraints_type & constraints);

      /**
        @brief Count of decompositions of @p mass within absolute tolerance @p error.

        Equivalent in result to @c getDecompositions(@p mass, @p error).size()
        but does not materialise the individual decompositions.

        When @c mass - error is non-positive, the integer-mass window
        starts at @c 1 rather than at a non-positive value.

        @param[in] mass  Target real-valued mass.
        @param[in] error Absolute mass tolerance.
        @return Number of decompositions whose reconstructed real
                parent mass is within @p error of @p mass.
      */
      number_of_decompositions_type getNumberOfDecompositions(double mass, double error);

private:
      /// Alphabet captured at construction time.
      Weights weights_;

      /// Minimum and maximum relative rounding errors of @ref weights_, captured at construction time.
      std::pair<double, double> rounding_errors_;

      /// Precision of @ref weights_, captured at construction time.
      double precision_;

      /// Lazy integer decomposer reused across calls; captured at construction time.
      std::shared_ptr<integer_decomposer_type> decomposer_;
    };

  } // namespace ims
} // namespace OpenMS

