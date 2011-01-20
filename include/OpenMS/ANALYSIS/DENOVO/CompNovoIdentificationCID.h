// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIDENTIFICATIONCID_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIDENTIFICATIONCID_H

// OpenMS includes
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h>

// stl includes
#include <vector>

namespace OpenMS
{
	/**
	  @brief  run with CompNovoIdentificationCID

		@htmlinclude OpenMS_CompNovoIdentificationCID.parameters	
	
		@ingroup Analysis_ID
	*/
	class OPENMS_DLLAPI CompNovoIdentificationCID : public CompNovoIdentificationBase
	{

		public:

			/** @name constructors and destructors
			 */
			//@{
			/// default constructor
			CompNovoIdentificationCID();
			
			/// copy constructor
			CompNovoIdentificationCID(const CompNovoIdentificationCID& source);
			
			/// destructor
			virtual ~CompNovoIdentificationCID();
			//@}
		
			///
			CompNovoIdentificationCID& operator = (const CompNovoIdentificationCID& source);

			/** @name Accessors
			 */
			//@{
			/// performs an ProteinIdentification run on a PeakMap
			void getIdentifications(std::vector<PeptideIdentification>& ids, const PeakMap& exp);

			/// performs an ProteinIdentification run on a PeakSpectrum
			void getIdentification(PeptideIdentification& id, const PeakSpectrum& CID_spec);
			//@}


			typedef CompNovoIonScoringBase::IsotopeType IsotopeType;
			typedef CompNovoIonScoringBase::IonScore IonScore;
			typedef CompNovoIdentificationBase::Permut Permut;
			
		protected:

			/// call the DAC algorithm for the subspectrum defined via left and right peaks and fill the set with candidates sequences
			void getDecompositionsDAC_(std::set<String>& sequences, Size left, Size right, DoubleReal peptide_weight, const PeakSpectrum& CID_orig_spec, Map<DoubleReal, IonScore>& CID_nodes);

			/// reduces the given number of permuts by scoring the perumtations to the CID and ETD spec
			void reducePermuts_(std::set<String>& permuts, const PeakSpectrum& CID_orig_spec, DoubleReal prefix, DoubleReal suffix);
			
			void updateMembers_();
			
			DoubleReal precursor_mass_tolerance_;
	};
}

#endif
