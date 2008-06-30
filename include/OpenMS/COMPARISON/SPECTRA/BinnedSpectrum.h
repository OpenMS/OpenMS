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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/SparseVector.h>
#include <OpenMS/CONCEPT/Exception.h>


#include <cmath>

namespace OpenMS
{

  /**
	@brief This is a binned representation of a PeakSpectrum

		@param 	sz the size of the bins and 
		@param	sp number of neighboring bins to both sides affected by a peak contribution
		@param	ps the peakspectrum, that shall be represented
		
		@p sz denotes the size of a bin in Th, thereby decideing the number of bins discretizeing 
		the spectrum. The peaks will be put in the respective bins and sum up inside. E.g. a BinnedSpectrum
		with binsize of 0.5 Th will have a peak at 100 Th in bin no. 200, a peak at 100.1 Th will be in bin 
		no. 201.
		@p sp denotes the number of neighboring bins to the left and the number of neighboring bins to the 
		right a peak is also added. 
		
		@ingroup SpectraComparison
  */
	
  class BinnedSpectrum : public MSSpectrum<>
  {
	  
  private:
	
	UInt binSpread_;
	Real binSize_;
	SparseVector<Real> bins_;
	
  public:

	/** 
		@brief 	Exception which is thrown if BinnedSpectrum bins are accessed and no PeakSpektrum has been 
				integrated yet i.e. bins_ is empty	
	*/
	class NoSpectrumIntegrated : public Exception::BaseException
	{
	public:
		NoSpectrumIntegrated(const char* file, int line, const char* function, const char* message ="BinnedSpectrum hasn't got a PeakSpectrum to base on yet") throw();

		virtual ~NoSpectrumIntegrated() throw();
	}; 
  
	typedef SparseVector<Real>::const_iterator const_bin_iterator;
	typedef SparseVector<Real>::iterator bin_iterator;
	
	/// default constructor
	BinnedSpectrum();
    
	/// constructor
	BinnedSpectrum(Real size, UInt spread, PeakSpectrum ps);
	
	/// copy constructor
	BinnedSpectrum(const BinnedSpectrum& source);

	/// destructor
	virtual ~BinnedSpectrum();
	
	/// assignment operator
	BinnedSpectrum& operator= (const BinnedSpectrum& source)
	{
		if (&source != this)
		{
			setBinSize(source.getBinSize());
			setBinSpread(source.getBinSpread());
			bins_ = source.getBins();
			MSSpectrum<>::operator=(source);
		}
		return *this;
	}	
	
	/// assignment operator for PeakSpectra
	BinnedSpectrum& operator= (const PeakSpectrum& source)
	{
		if (!MSSpectrum<>::operator==(source))
		{
			MSSpectrum<>::operator=(source);
			setBinning();
		}
		return *this;
	}	
	
	/// equality operator
	bool operator== (const BinnedSpectrum& rhs) const
	{
		return
			(MSSpectrum<>::operator==(rhs) &&
			rhs.getBinSize()==this->binSize_ &&
			rhs.getBinSpread()==this->binSpread_)
			;
	}
	
	/// inequality operator
	bool operator!= (const BinnedSpectrum& rhs) const
	{
		return !(operator==(rhs));
	}
		
	/// equality operator for PeakSpectra
	bool operator== (const PeakSpectrum& rhs) const
	{
		return MSSpectrum<>::operator==(rhs);
	}
	
	/// inequality operator for PeakSpectra
	bool operator!= (const PeakSpectrum& rhs) const
	{
		return !(operator==(rhs));
	}
    
	/// get the BinSize
	inline double getBinSize() const
	{
		return this->binSize_;
	}

	
	/// get the BinSpread
	inline UInt getBinSpread() const
	{
		return this->binSpread_;
	}
	
	/// get the BinNumber, number of Bins
	inline UInt getBinNumber() const
	{
		return (UInt)this->bins_.size();	
	}
	
	/// get the FilledBinNumber, number of filled Bins
	inline UInt getFilledBinNumber() const
	{
		return (UInt)this->bins_.nonzero_size();	
	}
	
	/// unmutable access to the Bincontainer
	inline const SparseVector<Real>& getBins() const throw (BinnedSpectrum::NoSpectrumIntegrated)
	{ 
		if(bins_.size() == 0)
		{
			throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return bins_; 
	}
	
	/// mutable access to the Bincontainer.
	inline SparseVector<Real>& getBins() throw (BinnedSpectrum::NoSpectrumIntegrated)
	{ 
		if(bins_.size() == 0)
		{
			try
			{
				this->setBinning();
			}
			catch(...)
			{
				throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}
		return bins_; 
	}
	
	/// returns the const begin iterator of the container
	inline const_bin_iterator begin() const 
	{ 
		return bins_.begin(); 
	}
	/// returns the const end iterator of the container
	inline const_bin_iterator end() const 
	{ 
		return bins_.end(); 
	}
	
	/// returns the begin iterator of the container
	inline bin_iterator begin() 
	{ 
		return bins_.begin(); 
	}
	/// returns the end iterator of the container
	inline bin_iterator end() 
	{ 
		return bins_.end(); 
	}
	
	/// sets the BinSize_ (and rebinnes)
	inline void setBinSize(double s) throw (BinnedSpectrum::NoSpectrumIntegrated)
	{
		if(this->binSize_ != s)
		{
			this->binSize_ = s;
			try
			{
				this->setBinning();
			}
			catch(...)
			{
				throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}		
		}
	}

	
	/// sets the BinSpread_ (and rebinnes)
	inline void setBinSpread(UInt s) throw (BinnedSpectrum::NoSpectrumIntegrated)
	{
		if(this->binSpread_ != s)
		{
			this->binSpread_ = s;
			try
			{
				this->setBinning();
			}
			catch(...)
			{
				throw BinnedSpectrum::NoSpectrumIntegrated(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}		
		}	
	}
	
	
	/// makes the binning: all Peaks of the containing PeakSpectrum are summed up in the bins corresponding to m/z ranges
	void setBinning() throw(BinnedSpectrum::NoSpectrumIntegrated);
	
	/// function to check compliance of two BinnedSpectrum objects
	bool checkCompliance(const BinnedSpectrum& bs) const;

	
  protected:
	// docu in base class
	virtual void clearChildIds_()
	{
		//TODO Persistence
	}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
