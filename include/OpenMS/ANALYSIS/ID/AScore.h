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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef  OPENMS_ANALYSIS_ID_ASCORE_H
#define  OPENMS_ANALYSIS_ID_ASCORE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
	class PeptideHit;
	/**
		@brief Implementation of the Ascore
			
		@todo Docu
	*/
	class OPENMS_DLLAPI AScore
	{
		public:
			///Default constructor
			AScore();
			///Destructor
			~AScore();

			/**
				@brief Computes the AScore. And returns a modified version of the input PeptideHit...

				@param	hit ....
				@param rp1d...
				@param fmt...fragment_mass_tolerance
			*/
			PeptideHit compute(PeptideHit& hit, RichPeakSpectrum& real_spectrum, DoubleReal fmt);
			
			///TODO
			Real computeCumulativeScore(UInt N,UInt n, Real p);
			
			///TODO
			void computeTwoHighestPeptides( std::vector< std::vector<Real> >& peptide_site_scores,Size& first, Size& second, Size& peak_depth);
			///TODO
			void compute_site_determining_ions(PeptideHit& hit, Size first,Size second, std::vector<RichPeakSpectrum>& site_determining_ions);
		private:
			/**
				@brief Computes number of matched ions between windows and the given spectrum...

				@param	th spectrum(should be theoretical)  ....
				@param windows  (should be 100 m/z)...
				@param fmt...fragment_mass_tolerance			
			*/
			Int numberOfMatchedIons_(const RichPeakSpectrum& th,const RichPeakSpectrum& windows ,Size depth, DoubleReal fmt);
		//	Int computeMatchedIons_();
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_ASCORE_H
