// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>


using namespace std;

namespace OpenMS
{
	BinnedSpectrum::BinnedSpectrum()
	   : MSSpectrum<>(), bin_spread_(1), bin_size_(2.0), bins_()
	{
	}

	BinnedSpectrum::BinnedSpectrum(Real size, UInt spread, PeakSpectrum ps)
	   : MSSpectrum<>(ps), bin_spread_(spread), bin_size_(size), bins_()
	{
		setBinning();
	}

	BinnedSpectrum::BinnedSpectrum(const BinnedSpectrum& source)
		: MSSpectrum<>(source), bin_spread_(source.getBinSpread()), bin_size_(source.getBinSize()), bins_(source.getBins())
	{
	}

	BinnedSpectrum::~BinnedSpectrum()
	{
	}


	//accessors and operators see .h file


	void BinnedSpectrum::setBinning()
	{
	 	if (this->MSSpectrum<>::empty())
	 	{
		 	throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		bins_.clear();

	 	//make all necessary bins accessible
		this->sortByPosition();
	 	bins_ = SparseVector<Real>((UInt)ceil(this->back().getMZ()/bin_size_) + bin_spread_ ,0,0);

		//put all peaks into bins
		UInt bin_number;
		for (Size i = 0; i < this->size(); ++i)
		{
			//bin_number counted form 0 -> floor
			bin_number = (UInt)floor(this->operator[](i).getMZ()/bin_size_);
			//e.g. bin_size_ = 1.5: first bin covers range [0,1.5] so peak at 1.5 falls in first bin (index 0)

			if(this->operator[](i).getMZ()/bin_size_ == (DoubleReal)bin_number)
			{
				--bin_number;
			}

			//add peak to corresponding bin
			bins_[bin_number] = bins_.at(bin_number) + this->operator[](i).getIntensity();

			//add peak to neighboring binspread many
			for (Size j = 0; j < bin_spread_; ++j)
			{
				bins_[bin_number+j+1] = bins_.at(bin_number+j+1) + this->operator[](i).getIntensity();
				// we are not in one of the first bins (0 to bin_spread)
				//not working:  if (bin_number-j-1 >= 0)
				if (bin_number >= j+1)
				{
					bins_[bin_number-j-1] = bins_.at(bin_number-j-1) + this->operator[](i).getIntensity();
				}
			}
		}

	}

	//yields false if given BinnedSpectrum size or spread differs from this one (comparing those might crash)
	bool BinnedSpectrum::checkCompliance(const BinnedSpectrum& bs) const
	{
		return
			(this->bin_size_ == bs.getBinSize()) &&
			(this->bin_spread_ == bs.getBinSpread());
	}

	BinnedSpectrum::NoSpectrumIntegrated::NoSpectrumIntegrated(const char* file, int line, const char* function, const char* message) throw()
          : BaseException(file, line, function, "BinnedSpectrum::NoSpectrumIntegrated", message)
	{
	}

	BinnedSpectrum::NoSpectrumIntegrated::~NoSpectrumIntegrated() throw()
	{
	}
}
