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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H
#define OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <string>
#include <vector>

namespace OpenMS
{
	/**
		 @brief Class for the enzymatic digestion of proteins
		 
		 @ingroup Chemistry
	*/
	class OPENMS_DLLAPI EnzymaticDigestion
	{
		public: 
			/// Possible enzymes for the digestion (adapt NamesOfEnzymes & getEnzymeByName() & nextCleavageSite_() if you add more enzymes here)
			enum Enzyme
			{
				TRYPSIN,
				SIZE_OF_ENZYMES
			};
			
			/// Names of the Enzymes
			static const std::string NamesOfEnzymes[SIZE_OF_ENZYMES];
			
			/// Default constructor
			EnzymaticDigestion();
			
			/// Returns the number of missed cleavages for the digestion
			Size getMissedCleavages() const;

			/// Sets the number of missed cleavages for the digestion (default is 0).
			void setMissedCleavages(Size missed_cleavages);	

			/// Returns the enzyme for the digestion
			Enzyme getEnzyme() const;
			
			/// Sets the enzyme for the digestion (default is TRYPSIN).
			void setEnzyme(Enzyme enzyme);		

			/// convert enzyme string name to enum
			/// returns SIZE_OF_ENZYMES if @p name is not valid
			Enzyme getEnzymeByName(const String& name);

			/// Performs the enzymatic digestion of a protein.
			void digest(const AASequence& protein, std::vector<AASequence>& output);
			
			/// Returns the number of peptides a digestion of @p protein would yield.
			Size peptideCount(const AASequence& protein);
			
		protected:
			/// Number of missed cleavages
			Size missed_cleavages_;
			/// Used enzyme
			Enzyme enzyme_;
			
			///moves the iterator @p it after the next cleavage site of the @p sequence
			void nextCleavageSite_(const AASequence& sequence, AASequence::ConstIterator& iterator);
	};

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ENZYMATICDIGESTION_H
