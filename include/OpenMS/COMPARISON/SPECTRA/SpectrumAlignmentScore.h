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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**
	  @brief Similarity score via spectra alignment

		This class implements a simple scoring based on the alignment of spectra. This alignment
		is implemented in the SpectrumAlignment class and performs a dynamic programming alignment
		of the peaks, minimizing the distances between the aligned peaks and maximizing the number
		of peak pairs. 

		The scoring is done via the simple formula score = sum / (sqrt(sum1 * sum2)). sum is the 
		product of the intensities of the aligned peaks, with the given exponent (default is 2). 
		sum1 and sum2 are the sum of the intensities squared for each peak of both spectra respectively.

		A binned version of this scoring is implemented in the ZhangSimilarityScoring class.
		 
		@htmlinclude OpenMS_SpectrumAlignmentScore.parameters
		
		@ingroup SpectraComparison
  */
	
  class OPENMS_DLLAPI SpectrumAlignmentScore : public PeakSpectrumCompareFunctor
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectrumAlignmentScore();

    /// copy constructor
    SpectrumAlignmentScore(const SpectrumAlignmentScore& source);

    /// destructor
    virtual ~SpectrumAlignmentScore();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SpectrumAlignmentScore& operator = (const SpectrumAlignmentScore& source);
	
		/// 
		double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;

		double operator () (const PeakSpectrum& spec) const;
		// @}

		// @name Accessors
		// @{
		///
    static PeakSpectrumCompareFunctor* create() { return new SpectrumAlignmentScore(); }

		///
		static const String getProductName()
		{
			return "SpectrumAlignmentScore";
		}
		// @}

		protected:
			
			/// returns the factor associated with the m/z tolerance and m/z difference of the peaks
			double getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian = false) const;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H
