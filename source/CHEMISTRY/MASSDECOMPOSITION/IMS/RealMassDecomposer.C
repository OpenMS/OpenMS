// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <iostream>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>

namespace OpenMS {

namespace ims {


RealMassDecomposer::RealMassDecomposer(const Weights& weights) :
  weights_(weights)
{

	rounding_errors_ = std::make_pair(weights.getMinRoundingError(), weights.getMaxRoundingError());
	precision_ = weights.getPrecision();
	decomposer_ = std::auto_ptr<integer_decomposer_type>(
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
         pos != decompositions.end();)
    {
			double parent_mass = weights_.getParentMass(*pos);
      if (fabs(parent_mass - mass) > error)
      {
				pos = decompositions.erase(pos);
			} else {
				++pos;
			}
		}
		all_decompositions_from_range.insert(all_decompositions_from_range.end(), 
                                         decompositions.begin(), decompositions.end());
	}
	return all_decompositions_from_range;
}

RealMassDecomposer::decompositions_type RealMassDecomposer::getDecompositions(double mass, double error,
                                      const constraints_type& constraints)
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
         pos != decompositions.end();) {
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
