// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <iostream>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>

namespace OpenMS
{
  namespace ims
  {

    RealMassDecomposer::RealMassDecomposer(const Weights & weights) :
      weights_(weights)
    {

      rounding_errors_ = std::make_pair(weights.getMinRoundingError(), weights.getMaxRoundingError());
      precision_ = weights.getPrecision();
      decomposer_ = std::shared_ptr<integer_decomposer_type>(
        new integer_decomposer_type(weights));
    }

    RealMassDecomposer::decompositions_type RealMassDecomposer::getDecompositions(double mass, double error)
    {
      // defines the range of integers to be decomposed
      integer_value_type start_integer_mass = static_cast<integer_value_type>(
        ceil((1 + rounding_errors_.first) * (mass - error) / precision_));
      integer_value_type end_integer_mass = static_cast<integer_value_type>(
        floor((1 + rounding_errors_.second) * (mass + error) / precision_));

      decompositions_type all_decompositions_from_range;

      // loops and finds decompositions for every integer mass,
      // then checks if real mass of decomposition lays in the allowed
      // error interval [mass-error; mass+error]
      for (integer_value_type integer_mass = start_integer_mass;
           integer_mass < end_integer_mass; ++integer_mass)
      {
        decompositions_type decompositions =
          decomposer_->getAllDecompositions(integer_mass);
        for (decompositions_type::iterator pos = decompositions.begin();
             pos != decompositions.end(); )
        {
          double parent_mass = weights_.getParentMass(*pos);
          if (fabs(parent_mass - mass) > error)
          {
            pos = decompositions.erase(pos);
          }
          else
          {
            ++pos;
          }
        }
        all_decompositions_from_range.insert(all_decompositions_from_range.end(),
                                             decompositions.begin(), decompositions.end());
      }
      return all_decompositions_from_range;
    }

    RealMassDecomposer::decompositions_type RealMassDecomposer::getDecompositions(double mass, double error,
                                                                                  const constraints_type & constraints)
    {

      // defines the range of integers to be decomposed
      integer_value_type start_integer_mass = static_cast<integer_value_type>(
        ceil((1 + rounding_errors_.first) * (mass - error) / precision_));
      integer_value_type end_integer_mass = static_cast<integer_value_type>(
        floor((1 + rounding_errors_.second) * (mass + error) / precision_));

      decompositions_type all_decompositions_from_range;

      // loops and finds decompositions for every integer mass,
      // then checks if real mass of decomposition lays in the allowed
      // error interval [mass-error; mass+error]
      for (integer_value_type integer_mass = start_integer_mass;
           integer_mass < end_integer_mass; ++integer_mass)
      {
        decompositions_type decompositions =
          decomposer_->getAllDecompositions(integer_mass);
        for (decompositions_type::iterator pos = decompositions.begin();
             pos != decompositions.end(); )
        {
          double parent_mass = weights_.getParentMass(*pos);
          if (fabs(parent_mass - mass) > error)
          {
            pos = decompositions.erase(pos);
          }
          else
          {
            bool to_erase = false;
            if (!constraints.empty())
            {
              for (constraints_type::const_iterator it =
                     constraints.begin(); it != constraints.end(); ++it)
              {
                if ((*pos)[it->first] < it->second.first ||
                    (*pos)[it->first] > it->second.second)
                {
                  to_erase = true;
                  break;
                }
              }
            }
            if (to_erase)
            {
              pos = decompositions.erase(pos);
            }
            else
            {
              ++pos;
            }
          }
        }
        all_decompositions_from_range.insert(all_decompositions_from_range.end(),
                                             decompositions.begin(), decompositions.end());
      }
      return all_decompositions_from_range;
    }

    RealMassDecomposer::number_of_decompositions_type RealMassDecomposer::getNumberOfDecompositions(double mass, double error)
    {
      // defines the range of integers to be decomposed
      integer_value_type start_integer_mass = static_cast<integer_value_type>(1);
      if (mass - error > 0)
      {
        start_integer_mass = static_cast<integer_value_type>(
          ceil((1 + rounding_errors_.first) * (mass - error) / precision_));
      }
      integer_value_type end_integer_mass = static_cast<integer_value_type>(
        floor((1 + rounding_errors_.second) * (mass + error) / precision_));

      number_of_decompositions_type number_of_decompositions = static_cast<number_of_decompositions_type>(0);

      // loops and finds decompositions for every integer mass,
      // then checks if real mass of decomposition lays in the allowed
      // error interval [mass-error; mass+error]
      for (integer_value_type integer_mass = start_integer_mass;
           integer_mass < end_integer_mass; ++integer_mass)
      {
        decompositions_type decompositions =
          decomposer_->getAllDecompositions(integer_mass);
        for (decompositions_type::iterator pos = decompositions.begin();
             pos != decompositions.end(); ++pos)
        {
          double parent_mass = weights_.getParentMass(*pos);
          if (fabs(parent_mass - mass) <= error)
          {
            ++number_of_decompositions;
          }
        }
      }
      return number_of_decompositions;
    }

  } // namespace ims
} // namespace OpenMS
