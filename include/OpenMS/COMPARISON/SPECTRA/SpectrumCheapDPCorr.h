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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{

  /**
  	@brief SpectrumCheapDPCorr calculates an optimal alignment on stick spectra
  
	  to save computing time only Peak Pairs that could get a score > 0 are considered
	  which Peak Pairs could get scores > 0 ? <br>
	  Peaks get a score depending on the difference in position and the heights of the peaks <br>
	  pairs with positions that differ more than some limit get score 0
  	
		@htmlinclude OpenMS_SpectrumCheapDPCorr.parameters

		@ingroup SpectraComparison
  */
	
  class OPENMS_DLLAPI SpectrumCheapDPCorr : public PeakSpectrumCompareFunctor
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectrumCheapDPCorr();

    /// copy constructor
    SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source);

    /// destructor
    virtual ~SpectrumCheapDPCorr();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SpectrumCheapDPCorr& operator = (const SpectrumCheapDPCorr& source);

		double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const;

		double operator () (const PeakSpectrum& a) const;
		// @}

		// @name Accessors
		// @{
		///
    static PeakSpectrumCompareFunctor* create() { return new SpectrumCheapDPCorr(); }

		///
		static const String getProductName()
		{
			return "SpectrumCheapDPCorr";
		}

    /// return consensus spectrum from last funtion call operator
    const PeakSpectrum& lastconsensus() const;

		///
		Map<UInt, UInt> getPeakMap() const;

    /// set weighting of the second spectrum for consensus from next function call operator
    void setFactor(double f);
		// @}

  private:

    /// O(n^2) dynamical programming
    double dynprog_(const PeakSpectrum&, const PeakSpectrum&, int, int, int, int) const;

    /// similarity of two peaks
    double comparepeaks_(double posa, double posb, double inta, double intb) const;

    static const String info_;

    /// consensus spectrum of the last comparison
    mutable PeakSpectrum lastconsensus_;

    /// should peaks with no alignment partner be kept in the consensus?
    bool keeppeaks_;

    /// weighting factor for the next consensus spectrum
    mutable double factor_;

		/// last peak map
		mutable Map<UInt, UInt> peak_map_;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H
