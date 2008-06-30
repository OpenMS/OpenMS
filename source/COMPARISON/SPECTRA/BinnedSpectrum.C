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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>


using namespace std;

namespace OpenMS
{  	
	BinnedSpectrum::BinnedSpectrum()
	   : MSSpectrum<>(), binSpread_(1), binSize_(2.0), bins_()
	{
	}
	
	BinnedSpectrum::BinnedSpectrum(Real size, UInt spread, PeakSpectrum ps)
	   : MSSpectrum<>(ps), binSpread_(spread), binSize_(size), bins_()
	{
		setBinning();
	}

	BinnedSpectrum::BinnedSpectrum(const BinnedSpectrum& source)
		: MSSpectrum<>(source), binSpread_(source.getBinSpread()), binSize_(source.getBinSize()), bins_(source.getBins())
	{
	}
	    
	BinnedSpectrum::~BinnedSpectrum()
	{
	}

	
	//accessors and operators see .h file
	
	
	void BinnedSpectrum::setBinning() throw (BinnedSpectrum::NoSpectrumIntegrated)
	{
	 	if (this->MSSpectrum<>::empty())
	 	{
		 	throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		} 
		this->bins_.clear();
	 	
	 	//make all necessary bins accessible
	 	bins_ = SparseVector<Real>((UInt)ceil(this->back().getMZ()/this->binSize_) + binSpread_ ,0,0);
	 	
		//put all peaks into bins 
		UInt bin_number;
		for (UInt i = 0; i < this->getContainer().size(); ++i)
		{
			//bin_number counted form 0 -> floor
			bin_number = (UInt)floor(this->getContainer()[i].getMZ()/this->binSize_);
			
			//add peak to corresponding bin
			bins_[bin_number] = bins_.at(bin_number) + this->getContainer()[i].getIntensity();
			
			//add peak to neighboring binspread many
			for (UInt j = 0; j < binSpread_; ++j)
			{
				bins_[bin_number+j+1] = bins_.at(bin_number+j+1) + this->getContainer()[i].getIntensity();
				// we are not in one of the first bins (0 to bin_spread)
				//not working:  if (bin_number-j-1 >= 0) 
				if (bin_number >= j+1) 
				{
					bins_[bin_number-j-1] = bins_.at(bin_number-j-1) + this->getContainer()[i].getIntensity();
				}
			}
		}
		
	}
	
	//yields false if given BinnedSpectrum size or spread differs from this one (comparing those might crash) 
	bool BinnedSpectrum::checkCompliance(const BinnedSpectrum& bs) const
	{
		return
			(this->binSize_ == bs.getBinSize()) &&
			(this->binSpread_ == bs.getBinSpread());
	}
			
	BinnedSpectrum::NoSpectrumIntegrated::NoSpectrumIntegrated(const char* file, int line, const char* function, const char* message) throw()
          : BaseException(file, line, function, "BinnedSpectrum::NoSpectrumIntegrated", message)
	{
	} 
	
	BinnedSpectrum::NoSpectrumIntegrated::~NoSpectrumIntegrated() throw()
	{
	}
}
