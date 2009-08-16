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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRASTSIMILARITYSCORE_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRASTSIMILARITYSCORE_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

namespace OpenMS
{

  /**
	  @brief Similarity score of SpectraST.
	  			 Unlike the other similarity scores this score is used for matching a spectrum against a whole library, although
	  			 the dot poduct seems to be an effective method for scoring on its own. For calculating the Spectrast score,
	  			 first preprocess the spectra if not already done. Transform them and calculate the dot product and the dot bias.
	  			 Afterwards get the best two hits and calculate delta_D. Now for every spectrum from the library you can calculate the final score.
	  			 

		The details of the score can be found in:
		H. Lam et al., Development and validation of a spectral library searching method for peptide identification from MS/MS,
		Proteomics, 7 , 655-667, 2007
		
		@ingroup SpectraComparison
  */
	
  class OPENMS_DLLAPI SpectraSTSimilarityScore : public PeakSpectrumCompareFunctor
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectraSTSimilarityScore();

    /// copy constructor
    SpectraSTSimilarityScore(const SpectraSTSimilarityScore& source);

    /// destructor
    virtual ~SpectraSTSimilarityScore();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SpectraSTSimilarityScore& operator = (const SpectraSTSimilarityScore& source);
	
		/**
			@brief: this function calculates only the dot product of the two spectrum like in SpectraST
		*/
		DoubleReal operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;
		
		DoubleReal operator() (const BinnedSpectrum& bin1,const BinnedSpectrum& bin2)	const;

		DoubleReal operator () (const PeakSpectrum& spec) const;
		
		/**
			@return true if spectrum passes filtering
		*/
		bool preprocess(PeakSpectrum& spec, Real remove_peak_intensity_threshold = 2.01, UInt cut_peaks_below = 1000, Size min_peak_number = 5, Size max_peak_number = 150);
		
		BinnedSpectrum transform(const PeakSpectrum& spec);
		/**
			@param dot_prdouct if -1 this value will be calculated as well.
		*/
		DoubleReal dot_bias(const BinnedSpectrum& bin1, const BinnedSpectrum& bin2, DoubleReal dot_product = -1) const;
		
		DoubleReal delta_D(DoubleReal top_hit, DoubleReal runner_up);
		
		DoubleReal compute_F(DoubleReal dot_product, DoubleReal delta_D, DoubleReal delta_bias);
		
		/**
			@brief: measures how much of the dot prudct is dominated by a few peaks
		*/
		// @}

		// @name Accessors
		// @{
		///
    static PeakSpectrumCompareFunctor* create() { return new SpectraSTSimilarityScore(); }

		///
		static const String getProductName()
		{
			return "SpectraSTSimilarityScore";
		}

		// @}

		protected:


  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRASTSIMILARTIYSCORE_H
