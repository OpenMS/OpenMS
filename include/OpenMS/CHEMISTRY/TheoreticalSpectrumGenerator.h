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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/CHEMISTRY/Residue.h>

namespace OpenMS
{
	class AASequence;
	
	/** 
		@brief Generates theoretical spectra with various options
		
		@todo replace SpectrumGenerator in CLUSTERING (Andreas)
		
		@ingroup Chemistry
	*/
	class TheoreticalSpectrumGenerator
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
			void getSpectrum(PeakSpectrum& spec, const AASequence& peptide, SignedInt charge = 1);

			/// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity
			void addPeaks(PeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, SignedInt charge = 1);

			/// adds the precursor peaks to the spectrum
			void addPrecursorPeaks(PeakSpectrum& spec, const AASequence& peptide, SignedInt charge = 1);

			/// mutable acces to the parameters
			Param& getParam();

			/// non-mutable access to the parameters
			const Param& getParam() const;

			/// set the parameters
			void setParam(const Param& param);
			//@}

		private:
			
			Param param_;

			Peak p_;
		};
}

#endif
