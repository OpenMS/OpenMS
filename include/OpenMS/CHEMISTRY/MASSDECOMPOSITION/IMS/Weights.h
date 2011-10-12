// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_WEIGHTS_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_WEIGHTS_H

#include <vector>
#include <ostream>

#include <OpenMS/config.h>

namespace OpenMS {

  namespace ims {

    /**
      @brief Represents a set of weights (double values and scaled with a certain
      precision their integer counterparts) with a quick access.

      Many algorithms can't work with real-valued alphabets and need integer
      weights. Those are usually obtained by dividing all alphabet masses by
      a "precision" parameter (0; 1) and rounding the result to integers.
      Class @c Weights allows access to the scaled masses (the weights) and to the
      original masses as well. Weights are cached and will not be recalculated on
      access.

      @note Words 'weights' and 'masses' were used deliberately for naming
      corresponding values to help users quickly understand how they can be used.
      Names were chosen due to similarity in 'weights' and 'masses' interconnection
      in this class with real life. In reality mass value is always the same, and
      weight value depends on circumstances. The same as here original double
      masses stay the same, and scaled integer weights depend on precision applied.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
    */
    class OPENMS_DLLAPI Weights
    {
    public:
      /**
        Type of integer values to be used.
      */
      typedef long unsigned int weight_type;

      /**
        Type of double values to be used.
      */
      typedef double alphabet_mass_type;

      /**
        Type of container to store integer values.
      */
      typedef std::vector<weight_type> weights_type;

      /**
        Type of container to store double values.
      */
      typedef std::vector<alphabet_mass_type> alphabet_masses_type;

      /**
        Type of container's size
      */
      typedef weights_type::size_type size_type;

      /**
        Empty constructor.
      */
      Weights() { }


      /**
        Constructor with double values and precision.

        @param masses Original double values to be scaled.
        @param precision Precision to scale double values.
      */
      Weights(const alphabet_masses_type& masses, alphabet_mass_type prec)
        : alphabet_masses_(masses),
          precision_(prec)
      {
        setPrecision(prec);
      }

      /**
        Copy constructor.

        @param other Weights to be copied.
      */
      Weights(const Weights& other) :
        alphabet_masses_(other.alphabet_masses_),
        precision_(other.precision_),
        weights_(other.weights_) { }

      /**
        Assignment operator.

        @param weights Weights to be assigned.
        @return Reference to this object.
      */
      Weights& operator =(const Weights& weights_);

      /**
        Gets size of a set of weights.

        @return Size of a set of weights.
      */
      size_type size() const { return weights_.size(); }

      /**
        Gets a scaled integer weight by index.

        @param i An index to access weights.
        @return An integer weight.
      */
      weight_type getWeight(size_type i) const { return weights_[i]; }

      /**
        Sets a new precision to scale double values to integer.
        *
        @param precision A new precision.
      */
      void setPrecision(alphabet_mass_type precision_);

      /**
        Gets precision.

        @return Precision to scale double values to integer.
      */
      alphabet_mass_type getPrecision() const { return precision_; }

      /**
        Operator to access weights by index.

        @param i An index to access weights.
        @return An integer weight.

        @see getWeight(size_type i)
      */
      weight_type operator [](size_type i) const { return weights_[i]; }

      /**
        Gets a last weight.

        @return a last weight.
      */
      weight_type back() const { return weights_.back(); }

      /**
        Gets an original (double) alphabet mass by index.

        @param i An index to access alphabet masses.
        @return A double alphabet mass.
      */
      alphabet_mass_type getAlphabetMass(size_type i) const
      { return alphabet_masses_[i]; }

      /**
        Returns a parent mass for a given \c decomposition
      */
      alphabet_mass_type getParentMass(const std::vector<unsigned int>& decomposition) const;

      /**
        Exchanges weight and mass at index1 with weight and mass at index2.

        @param index1 Index of weight and mass to be exchanged.
        @param index2 Index of weight and mass to be exchanged.
      */
      void swap(size_type index1, size_type index2);

      /**
        Divides the integer weights by their gcd. The precision is also
        adjusted.

        @return true if anything was changed, that is, if the gcd was &gt; 1.
        false if the gcd was already 1 or there are less than two weights.
      */
      bool divideByGCD();

      alphabet_mass_type getMinRoundingError() const;

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
      Prints weights to the stream @c os.

      @param os Output stream to which weights are written.
      @param weigths Weights to be written.
    */
    std::ostream& operator<<(std::ostream& os, const Weights& weights);

  } // namespace ims
} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_WEIGHTS_H
