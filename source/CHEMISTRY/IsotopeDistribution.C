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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Clemens Groepl, Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace std;

namespace OpenMS
{
	IsotopeDistribution::IsotopeDistribution()
		:	max_isotope_(0)
	{
		distribution_.push_back(make_pair<Size, double>(0, 1));
	}

	IsotopeDistribution::IsotopeDistribution(Size max_isotope)
		:	max_isotope_(max_isotope)
	{
		distribution_.push_back(make_pair<Size, double>(0, 1));
	}

	IsotopeDistribution::IsotopeDistribution(const IsotopeDistribution& isotope_distribution)
		:	max_isotope_(isotope_distribution.max_isotope_),
			distribution_(isotope_distribution.distribution_)
	{
	}

	IsotopeDistribution::~IsotopeDistribution()
	{
	}

	void IsotopeDistribution::setMaxIsotope(Size max_isotope)
	{
		max_isotope_ = max_isotope;
	}

	Size IsotopeDistribution::getMaxIsotope() const
	{
		return max_isotope_;
	}

	IsotopeDistribution& IsotopeDistribution::operator = (const IsotopeDistribution& iso)
	{
		if (this != &iso)
		{
			distribution_ = iso.distribution_;
			max_isotope_ = iso.max_isotope_;
		}
		return *this;
	}

	IsotopeDistribution IsotopeDistribution::operator + (const IsotopeDistribution& iso) const
	{
		ContainerType result;
		convolve_(result, distribution_, iso.distribution_);
		IsotopeDistribution result_iso;
		result_iso.setMaxIsotope(max_isotope_);
		result_iso.set(result);
		return result_iso;
	}

	IsotopeDistribution& IsotopeDistribution::operator += (const IsotopeDistribution& iso)
	{
		ContainerType result;
		convolve_(result, distribution_, iso.distribution_);
		distribution_ = result;
		return *this;
	}

	IsotopeDistribution& IsotopeDistribution::operator *= (Size factor)
	{
		ContainerType result;
		convolvePow_(result, distribution_, factor);
		distribution_ = result;
		return *this;
	}

	IsotopeDistribution IsotopeDistribution::operator * (Size factor) const
	{
		ContainerType result;
		convolvePow_(result, distribution_, factor);
		IsotopeDistribution result_iso;
		result_iso.setMaxIsotope(max_isotope_);
		result_iso.set(result);
		return result_iso;
	}

	void IsotopeDistribution::set(const ContainerType& distribution)
	{
		distribution_ = distribution;
	}

	const IsotopeDistribution::ContainerType& IsotopeDistribution::getContainer() const
	{
		return distribution_;
	}

	Size IsotopeDistribution::getMax() const
	{
		if (distribution_.size() == 0)
		{
			return 0;
		}
		return distribution_[distribution_.size()-1].first;
	}

	Size IsotopeDistribution::getMin() const
	{
		if (distribution_.size() == 0)
		{
			return 0;
		}
		return distribution_[0].first;
	}

	Size IsotopeDistribution::size() const
	{
		return distribution_.size();
	}

	void IsotopeDistribution::clear()
	{
		distribution_.clear();
		max_isotope_ = 0;
	}

	void IsotopeDistribution::estimateFromPeptideWeight(double average_weight)
	{
		const ElementDB* db = ElementDB::getInstance();

		vector<String> names;
		names.push_back("C");
		names.push_back("H");
		names.push_back("N");
		names.push_back("O");
		names.push_back("S");

		//Averagine element count divided by averagine weight
		vector<DoubleReal> factors;
		factors.push_back(4.9384/111.1254);
		factors.push_back(7.7583/111.1254);
		factors.push_back(1.3577/111.1254);
		factors.push_back(1.4773/111.1254);
		factors.push_back(0.0417/111.1254);

		//initialize distribution
		distribution_.clear();
		distribution_.push_back(make_pair(0u,1.0));

		for (Size i = 0; i != names.size(); ++i)
		{
			ContainerType single, conv_dist;
			//calculate distribution for single element
			ContainerType dist(db->getElement(names[i])->getIsotopeDistribution().getContainer());
			convolvePow_(single, dist, (Size)Math::round(average_weight * factors[i]));
			//convolve it with the existing distributions
			conv_dist = distribution_;
			convolve_(distribution_, single, conv_dist);
		}
	}

	bool IsotopeDistribution::operator == (const IsotopeDistribution& isotope_distribution) const
	{
		return 	max_isotope_ == isotope_distribution.max_isotope_ &&
						distribution_ == isotope_distribution.distribution_;
	}

	bool IsotopeDistribution::operator != (const IsotopeDistribution& isotope_distribution) const
	{
		return !(isotope_distribution == *this);
	}

	void IsotopeDistribution::convolve_(ContainerType& result, const ContainerType& left, const ContainerType& right) const
	{
		if (left.size() == 0 || right.size() == 0)
		{
			result.clear();
			return;
		}

		ContainerType::size_type r_max = left.size() + right.size() - 1;

		if ((ContainerType::size_type)max_isotope_ != 0 && r_max > (ContainerType::size_type)max_isotope_)
		{
			r_max = (ContainerType::size_type)max_isotope_;
		}

		result.resize(r_max);
    for (ContainerType::size_type i = 0; i != r_max; ++i)
    {
      result[i] = make_pair(left[0].first + right[0].first + i, 0);
    }

		// we loop backwards because then the small products tend to come first
		// (for better numerics)
		for (SignedSize i = left.size() - 1; i >= 0; --i)
		{
			for (SignedSize j = min<SignedSize>(r_max - i, right.size()) - 1; j >= 0; --j)
			{
				result[i+j].second += left[i].second * right[j].second;
			}
		}
	}

	void IsotopeDistribution::convolvePow_(ContainerType& result, const ContainerType& input, Size n) const
	{
		/*
		// my code
		ContainerType tmp, tmp_result;
		tmp.push_back(make_pair<Size, double>(0, 1));
		for (Size i=0; i!=n; ++i)
		{
			convolve(tmp_result, input, tmp);
			swap(tmp_result, tmp);
		}
		swap(tmp, result);
		*/


		// Clemens' code begin
		if (n == 1)
		{
			result = input;
			return;
		}

		Size log2n = 0;
		// modification by Chris to prevent infinite loop when n > 2^63
		if (n > (Size(1) << (std::numeric_limits<Size>::digits-1))) 
		{
			log2n = std::numeric_limits<Size>::digits;
		}
		else
		{
			// find binary logarithm of n
			for (; (Size(1) << log2n) < n; ++log2n) ;
		}
		

	  // get started
    if (n & 1)
		{
    	result = input;
		}
    else
		{
      result.clear();
      result.push_back(make_pair<Size, double>(0, 1.0));
    }

    ContainerType intermediate;

    // to avoid taking unneccessary squares, we check the loop condition
    // somewhere in the middle
    ContainerType convolution_power;
    convolveSquare_(convolution_power, input);
    for (Size i = 1; ; ++i)
    {
      if (n & (Size(1) << i))
      {
        convolve_(intermediate, result, convolution_power);
        swap(intermediate, result);
      }
      // check the loop condition
      if (i >= log2n) break;

      // prepare next round
      convolveSquare_(intermediate, convolution_power);
      swap(intermediate, convolution_power);
    }
		// Clemens' code end
	}

  void IsotopeDistribution::convolveSquare_(ContainerType& result, const ContainerType& input) const
  {
	  result.clear();
   	ContainerType::size_type r_max = 2 * input.size() - 1;

		if ((ContainerType::size_type)max_isotope_ != 0 && (ContainerType::size_type)(max_isotope_ + 1) < r_max)
		{
			r_max = (ContainerType::size_type)(max_isotope_ + 1);
		}

    result.resize(r_max);
		for (ContainerType::size_type i = 0; i != r_max; ++i)
		{
			result[i] = make_pair(2 * input[0].first + i, 0);
		}

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for (SignedSize i = input.size() - 1; i >= 0; --i)
		{
      for (SignedSize j = min<SignedSize>(r_max - i, input.size()) - 1; j >= 0; --j)
			{
        result[i+j].second += input[i].second * input[j].second;
      }
    }

    return;
  }


	void IsotopeDistribution::renormalize()
	{
		if (distribution_.size() != 0)
		{
			double sum(0);
			// loop backwards as most distributions contains a lot of small values at the end
			for (ConstIterator it = distribution_.end()-1; it != distribution_.begin(); --it)
			{
				sum += it->second;
			}
			sum += distribution_.begin()->second;

			for (Iterator it = distribution_.begin(); it != distribution_.end(); ++it)
			{
				it->second /= sum;
			}
		}
		return;
	}

	void IsotopeDistribution::trimRight(DoubleReal cutoff)
	{
		ContainerType::reverse_iterator riter = distribution_.rbegin();

		// loop from right to left until an entry is larger than the cutoff
		for ( ; riter != distribution_.rend(); riter++ )
		{
			if ( riter->second >= cutoff ) break;
		}
		// trim the container
		distribution_.resize( riter.base() - distribution_.begin() );
	}

	void IsotopeDistribution::trimLeft(DoubleReal cutoff)
	{
		for (ContainerType::iterator iter = distribution_.begin() ; iter != distribution_.end(); iter++ )
		{
			if ( iter->second >= cutoff )
			{
				distribution_.erase(distribution_.begin(),iter);
				break;
			}
		}
	}
}

