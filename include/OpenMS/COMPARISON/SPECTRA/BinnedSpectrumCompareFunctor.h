// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUMCOMPAREFUNCTOR_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUMCOMPAREFUNCTOR_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

#include <cmath>

namespace OpenMS
{

  /**
	
		@brief Base class for compare functors of BinnedSpectra
	
  		BinnedSpectrumCompareFunctor classes return a value for a pair of BinnedSpectrum objects (or a single one with itself).
  		Ideally the value should reflect the similarity of the pair. For methods of computing the similarity see the 
  		documetation of the concrete functors.
  		Functors normalized in the range [0,1] are identifiable at the set "normalized" parameter of the ParameterHandler
		
		@ingroup SpectraComparison
  */
  class OPENMS_DLLAPI BinnedSpectrumCompareFunctor : public DefaultParamHandler
  {
	  
  private:

  public:
    
    /** @brief Exception thrown if compared spectra are incompatible
	
		    the compared spectra have different settings in binsize and/or binspread 
		    due to which comparison would fail
    */
    class OPENMS_DLLAPI IncompatibleBinning
			: public Exception::BaseException
    {
    public:
      IncompatibleBinning(const char* file, int line, const char* function, const char* message
          = "compared spectra have different settings in binsize and/or binspread")  throw();
      virtual ~IncompatibleBinning() throw();
    };
	
    /// default constructor
    BinnedSpectrumCompareFunctor();

    /// copy constructor
    BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor& source);

    /// destructor
    virtual ~BinnedSpectrumCompareFunctor();

    /// assignment operator
    BinnedSpectrumCompareFunctor& operator = (const BinnedSpectrumCompareFunctor& source);

    /// function call operator, calculates the similarity of the given arguments
    virtual double operator () (const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const = 0;

	/// function call operator, calculates self similarity
	virtual double operator () (const BinnedSpectrum& spec) const = 0;

	/// registers all derived products 
	static void registerChildren();

	/// get the identifier for a DefaultParamHandler
	static const String getProductName()
	{
		return "BinnedSpectrumCompareFunctor";
	}
	
  };
  

}
#endif // OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUMCOMPAREFUNCTOR_H
