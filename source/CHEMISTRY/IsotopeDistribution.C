// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>

using namespace std;

namespace OpenMS
{
	IsotopeDistribution::IsotopeDistribution()
		:	max_isotope_(0)
	{
		distribution_.push_back(make_pair<UnsignedInt, double>(0, 1));
	}

	IsotopeDistribution::IsotopeDistribution(Size max_isotope)
		:	max_isotope_(max_isotope)
	{
		distribution_.push_back(make_pair<UnsignedInt, double>(0, 1));
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

	IsotopeDistribution::IsotopeDistribution& IsotopeDistribution::operator = (const IsotopeDistribution& iso)
	{
		if (this != &iso)
		{
			distribution_ = iso.distribution_;
			max_isotope_ = iso.max_isotope_;
		}
		return *this;
	}

	IsotopeDistribution::IsotopeDistribution IsotopeDistribution::operator + (const IsotopeDistribution& iso) const
	{
		ContainerType result;
		convolve_(result, distribution_, iso.distribution_);
		IsotopeDistribution result_iso;
		result_iso.setMaxIsotope(max_isotope_);
		result_iso.set(result);
		return result_iso;
	}

	IsotopeDistribution::IsotopeDistribution& IsotopeDistribution::operator += (const IsotopeDistribution& iso)
	{
		ContainerType result;
		convolve_(result, distribution_, iso.distribution_);
		distribution_ = result;
		return *this;
	}

	IsotopeDistribution::IsotopeDistribution& IsotopeDistribution::operator *= (UnsignedInt factor)
	{
		ContainerType result;
		convolvePow_(result, distribution_, factor);
		distribution_ = result;
		return *this;
	}

	IsotopeDistribution::IsotopeDistribution IsotopeDistribution::operator * (UnsignedInt factor) const
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

	UnsignedInt IsotopeDistribution::getMax() const
	{
		if (distribution_.size() == 0)
		{
			return 0;
		}
		return distribution_[distribution_.size()-1].first;
	}

	UnsignedInt IsotopeDistribution::getMin() const
	{
		if (distribution_.size() == 0)
		{
			return 0;
		}
		return distribution_[0].first;
	}

	UnsignedInt IsotopeDistribution::size() const
	{
		return distribution_.size();
	}

	void IsotopeDistribution::clear()
	{
		distribution_.clear();
		max_isotope_ = 0;
	}

	void IsotopeDistribution::estimateFromPeptideWeight(double weight)
	{
		// - there are 5.45 Carbons per Residue (average, assuming equal occuring frenquencies)
		// - about 53.6% of the monoisotopic weight of a peptide is from Carbon (average...)
		ContainerType C_dist;
		C_dist.push_back(make_pair<UnsignedInt, double>(12, 0.9893));
		C_dist.push_back(make_pair<UnsignedInt, double>(13, 0.0107));

		convolvePow_(distribution_, C_dist, Size((weight*0.464)/12.0));
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
		
		int r_max = min(left.size() + right.size() - 1, (ContainerType::size_type)max_isotope_);
		
		result.resize(r_max);
    for (int i = 0; i != r_max; ++i)
    {
      result[i] = make_pair<UnsignedInt, double>(left[0].first + right[0].first + i, 0);
    }

		// we loop backwards because then the small products tend to come first
		// (for better numerics)
		for (int i = left.size() - 1; i >= 0; --i)
		{
			for (int j = min(r_max - i, int(right.size())) - 1; j >= 0; --j)
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
		tmp.push_back(make_pair<UnsignedInt, double>(0, 1));
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

    // find binary logarithm of n
    Size log2n = 0;
    for (; (1U << log2n) < n; ++log2n);

	  // get started
    if (n & 1)
		{
    	result = input;
		}
    else 
		{
      result.clear();
      result.push_back(make_pair<UnsignedInt, double>(0, 1.0));
    }

    ContainerType intermediate;

    // to avoid taking unneccessary squares, we check the loop condition
    // somewhere in the middle
    ContainerType convolution_power;
    convolveSquare_(convolution_power, input);
    for (Size i = 1; ; ++i)
    {
      if (n & (1 << i))
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
   	int r_max = min(2 * input.size() - 1, (ContainerType::size_type)(max_isotope_ + 1));
    result.resize(r_max);
		for (int i = 0; i != r_max; ++i)
		{
			result[i] = make_pair<UnsignedInt, double>(2 * input[0].first + i, 0);
		}

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for (int i = input.size() - 1; i >= 0; --i) 
		{
      for (int j = min(r_max - i, int(input.size())) - 1; j >= 0; --j) 
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
	
	void IsotopeDistribution::trimRight()
	{
		ContainerType::reverse_iterator riter = distribution_.rbegin();
			
		// loop from right to left until an entry is larger than the cutoff
		for ( ; riter != distribution_.rend(); riter++ )
		{
			if ( riter->second >= trim_right_cutoff_ ) break;
		}
		// trim the container
		distribution_.resize ( riter.base() - distribution_.begin() );
	}	

}

