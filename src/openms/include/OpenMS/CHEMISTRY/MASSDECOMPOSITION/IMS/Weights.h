// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <vector>
#include <iosfwd>

#include <OpenMS/config.h>

namespace OpenMS
{
  namespace ims
  {
    /**
      @brief Pairs an alphabet of double-valued masses with the scaled
             integer weights derived from them.

      Many mass-decomposition algorithms only operate on integer weights;
      this class provides those integers alongside the original
      double-valued masses, derived once at construction (or on the next
      @ref setPrecision call) by computing
      @c round(mass / precision) per entry. Subsequent access is plain
      vector lookup -- no recomputation happens on read.

      Typical @c precision values are in @c (0, @c 1] (smaller values
      mean a finer-grained integer grid); the cpp itself does not enforce
      a range.

      @note The terms "weights" (the scaled integers) and "masses" (the
            original doubles) are intentionally kept distinct: the masses
            are invariants of the alphabet, while the weights depend on
            the chosen precision.

      @ingroup Analysis_DeNovo
    */
    class OPENMS_DLLAPI Weights
    {
public:
      /// Type of integer values to be used.
      typedef long unsigned int weight_type;

      /// Type of double values to be used.
      typedef double alphabet_mass_type;

      /// Type of container to store integer values.
      typedef std::vector<weight_type> weights_type;

      /// Type of container to store double values.
      typedef std::vector<alphabet_mass_type> alphabet_masses_type;

      /// Type of container's size
      typedef weights_type::size_type size_type;

      /**
        @brief Default constructor.

        Produces an empty instance; @ref size returns @c 0 until masses
        and a precision have been supplied.
      */
      Weights() {}

      /**
        @brief Construct from a set of masses and a precision.

        Stores @p masses as the alphabet masses and computes the integer
        weights via @c round(mass / precision) for each entry (so the
        @ref size of both the alphabet-mass and the weight container is
        @c masses.size() after the call).

        @param[in] masses    Original double-valued alphabet masses.
        @param[in] precision Precision used to scale the masses; smaller
                             values yield finer-grained integer weights.
      */
      Weights(const alphabet_masses_type & masses, alphabet_mass_type precision) :
        alphabet_masses_(masses),
        precision_(precision)
      {
        setPrecision(precision);
      }

      /**
        Copy constructor.

        @param[in] other Weights to be copied.
      */
      Weights(const Weights & other) :
        alphabet_masses_(other.alphabet_masses_),
        precision_(other.precision_),
        weights_(other.weights_) {}

      /**
        Assignment operator.

        @param[in] other Weights to be assigned.
        @return Reference to this object.
      */
      Weights & operator=(const Weights & other);

      /**
        Gets size of a set of weights.

        @return Size of a set of weights.
      */
      size_type size() const
      {
        return weights_.size();
      }

      /**
        Gets a scaled integer weight by index.

        @param[in] i An index to access weights.
        @return An integer weight.
      */
      weight_type getWeight(size_type i) const
      {
        return weights_[i];
      }

      /**
        @brief Replace the precision and recompute the integer weights from the stored alphabet masses.

        @param[in] precision New precision; the integer weights are
                             refreshed as @c round(mass / precision) for
                             every stored alphabet mass.
      */
      void setPrecision(alphabet_mass_type precision);

      /**
        Gets precision.

        @return Precision to scale double values to integer.
      */
      alphabet_mass_type getPrecision() const
      {
        return precision_;
      }

      /**
        Operator to access weights by index.

        @param[in] i An index to access weights.
        @return An integer weight.

        @see getWeight(size_type i)
      */
      weight_type operator[](size_type i) const
      {
        return weights_[i];
      }

      /**
        Gets a last weight.

        @return a last weight.
      */
      weight_type back() const
      {
        return weights_.back();
      }

      /**
        Gets an original (double) alphabet mass by index.

        @param[in] i An index to access alphabet masses.
        @return A double alphabet mass.
      */
      alphabet_mass_type getAlphabetMass(size_type i) const
      {
        return alphabet_masses_[i];
      }

      /**
        @brief Reconstruct a parent mass from a decomposition expressed as multiplicities of the alphabet.

        @param[in] decomposition Per-alphabet-entry multiplicities; must
                                 have the same length as the alphabet.
        @return @c sum(alphabet_mass[i] * decomposition[i]) using the
                stored original (double) masses.

        @throws Exception::InvalidParameter when @c decomposition.size()
                                            does not equal @ref size.
      */
      alphabet_mass_type getParentMass(const std::vector<unsigned int> & decomposition) const;

      /**
        @brief Exchange the alphabet mass @b and the integer weight at @p index1 with those at @p index2.

        Both the @c alphabet_mass entry and its corresponding integer
        @c weight at the two indices are swapped, so the mass / weight
        pairing is preserved.

        @param[in] index1 First index to swap.
        @param[in] index2 Second index to swap.
      */
      void swap(size_type index1, size_type index2);

      /**
        @brief Divide every integer weight by the greatest common divisor of all weights and scale the precision accordingly.

        For example, given alphabet weights @c 3.0, @c 5.0, @c 8.0 with
        precision @c 0.1, the integer weights would be @c 30, @c 50,
        @c 80. After calling this method the integer weights become
        @c 3, @c 5, @c 8 and the precision becomes @c 1.0 (since the
        gcd of @c 30, @c 50, @c 80 is @c 10).

        @return @c true when the integer weights and the precision were
                changed (gcd > 1); @c false when the gcd was already
                @c 1 or fewer than two weights are stored.
      */
      bool divideByGCD();

      /**
        @brief Largest negative relative rounding error introduced by the current precision.

        @return The minimum of @c (precision * weight[i] - mass[i]) / mass[i]
                over all entries with a negative value, or @c 0 if no
                entry's scaled weight under-estimates its mass.
      */
      alphabet_mass_type getMinRoundingError() const;

      /**
        @brief Largest positive relative rounding error introduced by the current precision.

        @return The maximum of @c (precision * weight[i] - mass[i]) / mass[i]
                over all entries with a positive value, or @c 0 if no
                entry's scaled weight over-estimates its mass.
      */
      alphabet_mass_type getMaxRoundingError() const;
private:
      /**
        Container to store original (double) alphabet masses.
      */
      alphabet_masses_type alphabet_masses_;

      /**
        Precision which is used to scale double values to integer.
      */
      alphabet_mass_type precision_;

      /**
        Container to store scaled integer weights.
      */
      weights_type weights_;
    };

    /**
      @brief Print every integer weight of @p weights to @p os, one per line.

      @param[out] os      Output stream.
      @param[in]  weights Weights container whose integer weights are written.
      @return Reference to @p os for chaining.
    */
    OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Weights & weights);

  } // namespace ims
} // namespace OpenMS

