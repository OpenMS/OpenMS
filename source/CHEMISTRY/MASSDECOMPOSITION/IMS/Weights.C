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

#include <cassert>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>

namespace OpenMS {

namespace ims {

/*
Weights::Weights(const alphabet_masses_type& masses, alphabet_mass_type precision) :
             alphabet_masses(masses),
             precision(precision) {
 setPrecision(precision);
}
*/

Weights& Weights::operator =(const Weights& other)
{
  if (this != &other)
  {
		alphabet_masses = other.alphabet_masses;
		precision = other.precision;
		weights = other.weights;
	}
	return *this;
}


void Weights::setPrecision(Weights::alphabet_mass_type precision)
{
	this->precision = precision;
	weights.clear();
	for (alphabet_masses_type::size_type i = 0;
       i < alphabet_masses.size(); ++i)
  {
		weights.push_back(static_cast<weight_type>(floor((alphabet_masses[i]
                                                     / precision) + 0.5)));
	}
}


void Weights::swap(size_type index1, size_type index2)
{
	weight_type weight = weights[index1];
	weights[index1] = weights[index2];
	weights[index2] = weight;

	alphabet_mass_type mass = alphabet_masses[index1];
	alphabet_masses[index1] = alphabet_masses[index2];
	alphabet_masses[index2] = mass;
}


Weights::alphabet_mass_type Weights::getParentMass(const std::vector<unsigned int>& decomposition) const
{
	alphabet_mass_type parent_mass = 0;
	
	assert(alphabet_masses.size() == decomposition.size());
  for (std::vector<unsigned int>::size_type i = 0; i < decomposition.size(); ++i)
  {
		parent_mass += alphabet_masses[i] * decomposition[i];
	}
	return parent_mass;	
}

/**
 * Divides the integer weights by their gcd. The precision is also adjusted.
 *
 * For example, given alphabet weights 3.0, 5.0, 8.0 with precision 0.1, the
 * integer weights would be 30, 50, 80. After calling this method, the new
 * weights are 3, 5, 8 with precision 1.0 (since the gcd of 30, 50, and 80
 * is 10).
 *
 * @return true if anything was changed, that is, if the gcd was &gt; 1.
 *         false if the gcd was already 1 or there are less than two weights.
*/
bool Weights::divideByGCD()
{
  if (weights.size() < 2)
  {
		return false;
	}
  weight_type d = Math::gcd(weights[0], weights[1]);
  for (weights_type::size_type i = 2; i < weights.size(); ++i)
  {
    d = Math::gcd(d, weights[i]);
		if (d == 1) {
			return false;
		}
	}
	// if we're here: d != 1

	precision *= d;

	// rescales the integer weights. Don't use setPrecision() here since
	// the result could be different due to rounding errors.
	for (weights_type::size_type i = 0; i < weights.size(); ++i) {
		weights[i] /= d;
	}
	return true;
}

Weights::alphabet_mass_type Weights::getMinRoundingError() const
{
	alphabet_mass_type min_error = 0;
	for (size_type i = 0; i < weights.size(); ++i) {
		alphabet_mass_type error = (precision * static_cast<alphabet_mass_type>(weights[i])
                                - alphabet_masses[i]) / alphabet_masses[i];
		if (error < 0 && error < min_error) {
			min_error = error;
		}
	}
	return min_error;
}

Weights::alphabet_mass_type Weights::getMaxRoundingError() const
{
	alphabet_mass_type max_error = 0;
	for (size_type i = 0; i < weights.size(); ++i) {
		alphabet_mass_type error = (precision * static_cast<alphabet_mass_type>(weights[i])
                                - alphabet_masses[i]) / alphabet_masses[i];
		if (error > 0 && error > max_error) {
			max_error = error;
		}
	}
	return max_error;
}

std::ostream& operator<<(std::ostream& os, const Weights& weights)
{
	for (Weights::size_type i = 0; i < weights.size(); ++i ) {
		os << weights.getWeight(i) << std::endl;
	}
	return os;
}

} // namespace ims

} // namespace OpenMS
