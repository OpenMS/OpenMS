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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARITYSCORE_H
#define OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARITYSCORE_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**
	  @brief Similarity score of Zhang

		The details of the score can be found in:
		Z. Zhang, Prediction of Low-Energy Collision-Induced Dissociation Spectra of Peptides,
		Anal. Chem., 76 (14), 3908 - 3922, 2004

		@htmlinclude OpenMS_ZhangSimilarityScore.parameters
		
		@ingroup SpectraComparison
  */
	
  class OPENMS_DLLAPI ZhangSimilarityScore : public PeakSpectrumCompareFunctor
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    ZhangSimilarityScore();

    /// copy constructor
    ZhangSimilarityScore(const ZhangSimilarityScore& source);

    /// destructor
    virtual ~ZhangSimilarityScore();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    ZhangSimilarityScore& operator = (const ZhangSimilarityScore& source);
	
		/// 
		double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;

		double operator () (const PeakSpectrum& spec) const;
		// @}

		// @name Accessors
		// @{
		///
    static PeakSpectrumCompareFunctor* create() { return new ZhangSimilarityScore(); }

		///
		static const String getProductName()
		{
			return "ZhangSimilarityScore";
		}

		// @}

		protected:

      /// returns the factor associated with the m/z tolerance and m/z difference of the peaks
      double getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian = false) const;


  };

}
#endif //OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARTIYSCORE_H
