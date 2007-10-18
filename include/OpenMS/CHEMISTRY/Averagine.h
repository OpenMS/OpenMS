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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_AVERAGINE_H
#define OPENMS_CHEMISTRY_AVERAGINE_H


#include <OpenMS/KERNEL/RawDataPoint1D.h>
#include <vector>

#ifndef GAY_CONSTANTS
#define GAY_CONSTANTS 1

#define M01 0.0576
#define M02 -0.2553
#define M03 0.5827
#define M04 -0.8436
#define M05 0.6182

#define M11 -0.1207
#define M12 0.4688
#define M13 -0.7373
#define M14 0.3845
#define M15 0.2701

#define M21 0.0569
#define M22 -0.1173
#define M23 -0.0838
#define M24 0.3068
#define M25 0.0865

#define M31 0.0296
#define M32 -0.1293
#define M33 0.1330
#define M34 0.1151
#define M35 0.0204

#define M41 -0.0101
#define M42 -0.0118
#define M43 0.0798
#define M44 0.0292
#define M45 0.0040

#define M51 -0.0132
#define M52 0.0448
#define M53 0.0256
#define M54 0.0079
#define M55 0.0008

#endif


namespace OpenMS
{

	/** 
		@brief A class modelling the isotopic distribution in dependence of the mass m. 
		
		The implementation follows the ideas and techniques described in:
		S. Gay et al.: Modeling peptide mass fingerprinting data using the atomic composition of peptides.
		Electrophoresis (20), 3527-34, 1999.
		
		@todo Merge with IsotopeDistribution? (Rene, Marc, Clemens, Andreas)
	*/
	class Averagine
	{
		public:
	
			/** 
				@brief Returns the peak distribution for a given mass. 
	 			
	 			@param mz_pos The monoistopic position for which an isotopic pattern should be derived.
	 			@param charge The corresponding charge the model should implement.
	 			@param pattern_extend An optional STL pair of doubles. If provided, the model returns
	 			  as first value the estimated m/z position at which the pattern appears.
	 			  Consequently, the second value resembles the "end position" of the pattern.
	 			@param cut_off Determines when a pattern should be considered as vanished. Hence, a value of 0.01 indicates
	 			  a cutoff a soon as the signal's intensity drops below 1% of the sum over all present peak intensities. 
	 		*/	
			static std::vector<RawDataPoint1D> getModel (const double mz_pos, const unsigned int charge, std::pair<double, double>* pattern_extend=NULL, const double cut_off=0.01) throw ();
	
		protected:
	
			/** Internally used constants according to Gay et al. */
			static double gay_constants[6][5]; 
	};

} // namespace OpenMS

#endif

