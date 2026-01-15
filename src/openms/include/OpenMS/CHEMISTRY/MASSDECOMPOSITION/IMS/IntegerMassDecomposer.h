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
#include <utility>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/MassDecomposer.h>

#include <OpenMS/MATH/MathFunctions.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief Implements @c MassDecomposer interface using algorithm and data
      structures described in paper "Efficient Mass Decomposition"
      S. Böcker, Z. Lipták, ACM SAC-BIO, 2004 doi:10.1145/1066677.1066715.

      The main idea is instead of using the classical dynamic programming
      algorithm, store the residues of the smallest decomposable numbers
      for every modulo of the smallest alphabet mass.

      @tparam ValueType Type of values to be decomposed.
      @tparam DecompositionValueType Type of decomposition elements.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
      @author Marcel Martin <Marcel.Martin@CeBiTec.Uni-Bielefeld.DE>
      @author Henner Sudek <Henner.Sudek@CeBiTec.Uni-Bielefeld.DE>
    */
    template <typename ValueType = long unsigned int,
              typename DecompositionValueType = unsigned int>
    class IntegerMassDecomposer :
      public MassDecomposer<ValueType, DecompositionValueType>
    {
public:
      /// Type of value to be decomposed.
      typedef typename MassDecomposer<ValueType, DecompositionValueType>::value_type value_type;

      /// Type of decomposition value.
      typedef typename MassDecomposer<ValueType, DecompositionValueType>::decomposition_value_type decomposition_value_type;

      /// Type of decomposition.
      typedef typename MassDecomposer<ValueType, DecompositionValueType>::decomposition_type decomposition_type;

      /// Type of container for many decompositions.
      typedef typename MassDecomposer<ValueType, DecompositionValueType>::decompositions_type decompositions_type;

      /// Type of decomposition's size.
      typedef typename decomposition_type::size_type size_type;

      /**
        Constructor with weights.

        @param alphabet Weights over which masses to be decomposed.
      */
      explicit IntegerMassDecomposer(const Weights & alphabet);

      /**
        Returns true if decomposition over the @c mass exists, otherwise - false.

        @param mass Mass to be decomposed.
        @return true if decomposition over a given mass exists, otherwise - false.
      */
      bool exist(value_type mass) override;

      /**
        Gets one possible decomposition for @c mass.

        @param mass Mass to be decomposed.
        @return One possible decomposition for a given mass.
      */
      decomposition_type getDecomposition(value_type mass) override;

      /**
        Gets all possible decompositions for @c mass.

        @param mass Mass to be decomposed.
        @return All possible decompositions for a given mass.
      */
      decompositions_type getAllDecompositions(value_type mass) override;

      /**
        Gets number of all possible decompositions for a given @p mass.
        Since using getAllDecomposition() the usage of this function could
        be @b consuming.

        @param mass Mass to be decomposed
        @return number of decompositions for a given mass.
      */
      decomposition_value_type getNumberOfDecompositions(value_type mass) override;

private:

      /**
        Type of witness vector.
      */
      typedef std::vector<std::pair<size_type, decomposition_value_type> > witness_vector_type;

      /**
        Type of rows of residues table.
      */
      typedef std::vector<value_type> residues_table_row_type;

      /**
        Type of the residues table.
      */
      typedef std::vector<residues_table_row_type> residues_table_type;

      /**
        Weights over which the mass is to be decomposed.
      */
      Weights alphabet_;

      /**
        Table with the residues of the smallest decomposable numbers over
        every modulo of the smallest alphabet mass are stored.
        Corresponds to the Extended Residue Table in the paper.
      */
      residues_table_type ertable_;

      /**
        List of the least common multiples. Corresponds to the lcm data structure
        in the paper.
      */
      residues_table_row_type lcms_;

      /**
        List of the counters for the least common multiples that store the number
        how often the smallest alphabet mass fits into the corresponding least
        common multiple(lcm).
      */
      residues_table_row_type mass_in_lcms_;

      /**
        Represents an infinite value.
      */
      value_type infty_;

      /**
        List of the witnesses that is used to find one mass decomposition.
        Corresponds to the witness vector w in the paper.
      */
      witness_vector_type witness_vector_;

      /**
        Fills the extended residues table.
      */
      void fillExtendedResidueTable_(const Weights & _alphabet, residues_table_row_type & _lcms,
                                     residues_table_row_type & _mass_in_lcms, const value_type _infty,
                                     witness_vector_type & _witness_vector, residues_table_type & _ertable);

      /**
        Collects decompositions for @c mass by recursion.

        @param mass Mass to be decomposed.
        @param alphabetMassIndex An index of the mass in alphabet that is used on this step of recursion.
        @param decomposition Decomposition which is calculated on this step of recursion.
        @param decompositionsStore Container where decompositions are collected.
      */
      void collectDecompositionsRecursively_(value_type mass, size_type alphabetMassIndex,
                                             decomposition_type decomposition, decompositions_type & decompositionsStore);
    };


    template <typename ValueType, typename DecompositionValueType>
    IntegerMassDecomposer<ValueType, DecompositionValueType>::IntegerMassDecomposer(
      const Weights & alphabet) :
      alphabet_(alphabet)
    {

      lcms_.resize(alphabet.size());
      mass_in_lcms_.resize(alphabet.size());

      infty_ = alphabet.getWeight(0) * alphabet.getWeight(alphabet.size() - 1);

      fillExtendedResidueTable_(alphabet, lcms_, mass_in_lcms_, infty_, witness_vector_, ertable_);

    }

    template <typename ValueType, typename DecompositionValueType>
    void IntegerMassDecomposer<ValueType, DecompositionValueType>::fillExtendedResidueTable_(
      const Weights & _alphabet, residues_table_row_type & _lcms, residues_table_row_type & _mass_in_lcms,
      const value_type _infty, witness_vector_type & _witnessVector, residues_table_type & _ertable)
    {

      if (_alphabet.size() < 2)
      {
        return;
      }
      // caches the most often used mass - smallest mass
      value_type smallestMass = _alphabet.getWeight(0), secondMass = _alphabet.getWeight(1);

      // initializes table: infinity everywhere except in the first field of every column
      _ertable.reserve(_alphabet.size());
      _ertable.assign(_alphabet.size(), std::vector<value_type>(smallestMass, _infty));

      for (size_type i = 0; i < _alphabet.size(); ++i)
      {
        _ertable[i][0] = 0;
      }

      // initializes witness vector
      _witnessVector.resize(smallestMass);

      // fills second column (the first one is already correct)
      size_type it_inc = secondMass % smallestMass, witness = 1;
      //typename residues_table_row_type::iterator it = _ertable[1].begin() + it_inc;
      value_type mass = secondMass;
      // initializes counter to create a witness vector
      decomposition_value_type counter = 0;
      size_type it_i = it_inc;
      while (it_i != 0)
      {
        _ertable[1][it_i] = mass;
        mass += secondMass;
        ++counter;
        _witnessVector[it_i] = std::make_pair(witness, counter);
        //std::cerr << "BLA: " << counter << " " << &_ertable[1][0] << " " << it - _ertable[1].begin() << " " << _ertable[1].size() << std::endl;
        it_i += it_inc;
        if (it_i >= _ertable[1].size())
        {
          it_i -= _ertable[1].size();
        }
      }
      // fills cache variables for i==1
      value_type tmp_d = Math::gcd(smallestMass, secondMass);
      _lcms[1] = secondMass * smallestMass / tmp_d;
      _mass_in_lcms[1] = smallestMass / tmp_d;

      // fills remaining table. i is the column index.
      for (size_type i = 2; i < _alphabet.size(); ++i)
      {
        // caches often used i-th alphabet mass
        value_type currentMass = _alphabet.getWeight(i);

        value_type d = Math::gcd(smallestMass, currentMass);

        // fills cache for various variables.
        // note that values for i==0 are never assigned since they're unused anyway.
        _lcms[i] = currentMass * smallestMass / d;
        _mass_in_lcms[i] = smallestMass / d;

        // Nijenhuis' improvement: Is currentMass composable with smaller alphabet?
        if (currentMass >= _ertable[i - 1][currentMass % smallestMass])
        {
          _ertable[i] = _ertable[i - 1];
          continue;
        }

        const residues_table_row_type & prev_column = _ertable[i - 1];
        residues_table_row_type & cur_column = _ertable[i];

        if (d == 1)
        {
          // This loop is for the case that the gcd is 1. The optimization used below
          // is not applicable here.

          // p_inc is used to change residue (p) efficiently
          size_type p_inc = currentMass % smallestMass;

          // n is the value that will be written into the table
          value_type n = 0;
          // current residue (in paper variable 'r' is used)
          size_type p = 0;
          // counter for creation of witness vector
          decomposition_value_type local_counter = 0;

          for (size_type m = smallestMass; m > 0; --m)
          {
            n += currentMass;
            p += p_inc;
            ++local_counter;
            if (p >= smallestMass)
            {
              p -= smallestMass;
            }
            if (n > prev_column[p])
            {
              n = prev_column[p];
              local_counter = 0;
            }
            else
            {
              _witnessVector[p] = std::make_pair(i, local_counter);
            }
            cur_column[p] = n;
          }
        }
        else
        {
          // If we're here, the gcd is not 1. We can use the following cache-optimized
          // version of the algorithm. The trick is to put the iteration over all
          // residue classes into the _inner_ loop.
          //
          // One could see it as going through one column in blocks which are gcd entries long.
          size_type cur = currentMass % smallestMass;
          size_type prev = 0;
          size_type p_inc = cur - d;
          // counters for creation of one witness vector
          std::vector<decomposition_value_type> counters(smallestMass);

          // copies first block from prev_column to cur_column
          for (size_type j = 1; j < d; ++j)
          {
            cur_column[j] = prev_column[j];
          }

          // first loop: goes through all blocks, updating cur_column for the first time.
          for (size_type m = smallestMass / d; m > 1; m--)
          {
            // r: current residue class
            for (size_type r = 0; r < d; r++)
            {

              ++counters[cur];
              if (cur_column[prev] + currentMass > prev_column[cur])
              {
                cur_column[cur] = prev_column[cur];
                counters[cur] = 0;
              }
              else
              {
                cur_column[cur] = cur_column[prev] + currentMass;
                _witnessVector[cur] = std::make_pair(i, counters[cur]);
              }

              prev++;
              cur++;
            }

            prev = cur - d;

            // this does: cur = (cur + currentMass) % smallestMass - d;
            cur += p_inc;
            if (cur >= smallestMass)
            {
              cur -= smallestMass;
            }
          }

          // second loop:
          bool cont = true;
          while (cont)
          {
            cont = false;
            prev++;
            cur++;
            ++counters[cur];
            for (size_type r = 1; r < d; ++r)
            {
              if (cur_column[prev] + currentMass < cur_column[cur])
              {
                cur_column[cur] = cur_column[prev] + currentMass;
                cont = true;
                _witnessVector[cur] = std::make_pair(i, counters[cur]);
              }
              else
              {
                counters[cur] = 0;
              }
              prev++;
              cur++;
            }

            prev = cur - d;

            cur += p_inc;
            if (cur >= smallestMass)
            {
              cur -= smallestMass;
            }
          }
        }

      }
    }

    template <typename ValueType, typename DecompositionValueType>
    bool IntegerMassDecomposer<ValueType, DecompositionValueType>::
    exist(value_type mass)
    {

      value_type residue = ertable_.back().at(mass % alphabet_.getWeight(0));
      return residue != infty_ && mass >= residue;
    }

    template <typename ValueType, typename DecompositionValueType>
    typename IntegerMassDecomposer<ValueType, DecompositionValueType>::decomposition_type
    IntegerMassDecomposer<ValueType, DecompositionValueType>::getDecomposition(value_type mass)
    {

      decomposition_type decomposition;
      if (!this->exist(mass))
      {
        return decomposition;
      }

      decomposition.reserve(alphabet_.size());
      decomposition.resize(alphabet_.size());

      // initial mass residue: in FIND-ONE algorithm in paper corresponds variable "r"
      value_type r = mass % alphabet_.getWeight(0);
      value_type m = ertable_.back().at(r);

      decomposition.at(0) = static_cast<decomposition_value_type>
                            ((mass - m) / alphabet_.getWeight(0));

      while (m != 0)
      {
        size_type i = witness_vector_.at(r).first;
        decomposition_value_type j = witness_vector_.at(r).second;
        decomposition.at(i) += j;
        if (m < j * alphabet_.getWeight(i))
        {
          break;
        }
        m -= j * alphabet_.getWeight(i);
        r = m % alphabet_.getWeight(0);
      }
      return decomposition;
    }

    template <typename ValueType, typename DecompositionValueType>
    typename IntegerMassDecomposer<ValueType, DecompositionValueType>::decompositions_type
    IntegerMassDecomposer<ValueType, DecompositionValueType>::getAllDecompositions(value_type mass)
    {
      decompositions_type decompositionsStore;
      decomposition_type decomposition(alphabet_.size());
      collectDecompositionsRecursively_(mass, alphabet_.size() - 1, decomposition, decompositionsStore);
      return decompositionsStore;
    }

    template <typename ValueType, typename DecompositionValueType>
    void IntegerMassDecomposer<ValueType, DecompositionValueType>::
    collectDecompositionsRecursively_(value_type mass, size_type alphabetMassIndex,
                                      decomposition_type decomposition, decompositions_type & decompositionsStore)
    {
      if (alphabetMassIndex == 0)
      {
        value_type numberOfMasses0 = mass / alphabet_.getWeight(0);
        if (numberOfMasses0 * alphabet_.getWeight(0) == mass)
        {
          decomposition[0] = static_cast<decomposition_value_type>(numberOfMasses0);
          decompositionsStore.push_back(decomposition);
        }
        return;
      }

      // tested: caching these values gives us 15% better performance, at least
      // with aminoacid-mono.masses
      const value_type lcm = lcms_[alphabetMassIndex];
      const value_type mass_in_lcm = mass_in_lcms_[alphabetMassIndex]; // this is alphabet mass divided by gcd

      value_type mass_mod_alphabet0 = mass % alphabet_.getWeight(0); // trying to avoid modulo
      const value_type mass_mod_decrement = alphabet_.getWeight(alphabetMassIndex) % alphabet_.getWeight(0);

      for (value_type i = 0; i < mass_in_lcm; ++i)
      {
        // here is the conversion from value_type to decomposition_value_type
        decomposition[alphabetMassIndex] = static_cast<decomposition_value_type>(i);

        // this check is needed because mass could have unsigned type and after reduction on i*alphabetMass will be still be positive but huge
        // and that will end up in infinite loop
        if (mass < i * alphabet_.getWeight(alphabetMassIndex))
        {
          break;
        }

        // r: current residue class. will stay the same in the following loop
        value_type r = ertable_[alphabetMassIndex - 1][mass_mod_alphabet0];

        // TODO: if infty was std::numeric_limits<...>... the following 'if' would not be necessary
        if (r != infty_)
        {
          for (value_type m = mass - i * alphabet_.getWeight(alphabetMassIndex); m >= r; m -= lcm)
          {
            // the condition of the 'for' loop (m >= r) and decrementing the mass
            // in steps of the lcm ensures that m is decomposable. Therefore
            // the recursion will result in at least one witness.
            collectDecompositionsRecursively_(m, alphabetMassIndex - 1, decomposition, decompositionsStore);
            decomposition[alphabetMassIndex] += mass_in_lcm;
            // this check is needed because mass could have unsigned type and after reduction on i*alphabetMass will be still be positive but huge
            // and that will end up in infinite loop
            if (m < lcm)
            {
              break;
            }
          }
        }
        // subtle way of changing the modulo, instead of plain calculation it from (mass - i*currentAlphabetMass) % alphabetMass0 every time
        if (mass_mod_alphabet0 < mass_mod_decrement)
        {
          mass_mod_alphabet0 += alphabet_.getWeight(0) - mass_mod_decrement;
        }
        else
        {
          mass_mod_alphabet0 -= mass_mod_decrement;
        }
      }

    }

    /**
      Gets number of all possible decompositions for a given @p mass.
      Since using getAllDecomposition() the usage of this function could
      be @b consuming.

      Needs the @p mass to be decomposed and returns the number of decompositions for that mass.
    */
    template <typename ValueType, typename DecompositionValueType>
    typename IntegerMassDecomposer<ValueType, DecompositionValueType>::decomposition_value_type IntegerMassDecomposer<ValueType,
                                                                                                                      DecompositionValueType>::getNumberOfDecompositions(value_type mass)
    {
      return static_cast<typename IntegerMassDecomposer<ValueType, DecompositionValueType>::decomposition_value_type>(getAllDecompositions(mass).size());
    }

  } // namespace ims
} // namespace OpenMS

