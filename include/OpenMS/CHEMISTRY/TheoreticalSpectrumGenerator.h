// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
	class AASequence;
	
	/** 
		@brief Generates theoretical spectra with various options
   	
   	@htmlinclude OpenMS_TheoreticalSpectrumGenerator.parameters
		
		@ingroup Chemistry
	*/
	class TheoreticalSpectrumGenerator : public DefaultParamHandler
	{
		public:
		
			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			TheoreticalSpectrumGenerator();

			/// copy constructor
			TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source);

			/// destructor
			virtual ~TheoreticalSpectrumGenerator();
			//@}
			
			/// assignment operator 
			TheoreticalSpectrumGenerator& operator = (const TheoreticalSpectrumGenerator& tsg);

			/** @name Acessors
			 */
			//@{
			/// returns a spectrum with b and y peaks
			void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1);

			/// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity
			void addPeaks(RichPeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, Int charge = 1);

			/// adds the precursor peaks to the spectrum
			void addPrecursorPeaks(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1);
			//@}

		private:
			
			RichPeak1D p_;
		};
}

#endif
