// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2004 -- Oliver Kohlbacher, Knut Reinert
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
// $Id:$
// $Author:$
// $Maintainer:$
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/CHEMISTRY/PeptideSequence.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/** \brief Generates theoretical spectra with various options
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
			TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator&);

			/// destructor
			virtual ~TheoreticalSpectrumGenerator();
			//@}
			
			/// assignment operator 
			TheoreticalSpectrumGenerator operator = (const TheoreticalSpectrumGenerator& tsg);

			/// returns a spectrum with b and y peaks
			PeakSpectrum getSpectrum(const PeptideSequence& peptide);

			/// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity
			void addPeaks(PeakSpectrum& spectrum, const PeptideSequence& peptide,
														Residue::ResidueType res_type, SignedInt charge = 1, double intensity = 1);

			/// returns a vector with sequences of the ions
			std::vector<PeptideSequence> getIons() const;

			/// returns true if loss peak adding is enabled 
			bool getAddLosses() const;
			
			/// setter to toggle loss peak adding
			void setAddLosses(bool add_losses);
		
			/// returns true if isotope peaks adding is enabled
			bool getAddIsotopes() const;

			/// setter to toggle isotope peaks adding
			void setAddIsotopes(bool add_isotopes);

			/// return true if meta info support is enabled
			bool getAddMetaInfo() const;

			/// setter to toggle meta info support
			void setAddMetaInfo(bool meta_info);

			/// returns the max isotope which is added via isotope peak support (if enabled)
			UnsignedInt getMaxIsotope() const;

			/// set the max peak of the isotope support
			void setMaxIsotope(UnsignedInt max_isotope);

		private:
			
			bool add_losses_;

			bool add_isotopes_;

			bool add_metainfo_;

			UnsignedInt max_isotope_;

			std::vector<PeptideSequence> ions_;
	};
}

#endif
