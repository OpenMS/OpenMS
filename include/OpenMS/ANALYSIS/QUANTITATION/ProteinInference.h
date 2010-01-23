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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_PROTEININFERENCE_H
#define OPENMS_ANALYSIS_QUANTITATION_PROTEININFERENCE_H

#include <vector>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{
	/**
		@brief [experimental class] given a peptide quantitation, infer corresponding protein quantities
		
		Infers protein ratios from peptide ratios (currently using unique peptides only).
		Use the IDMapper class to add protein and peptide information to a 
		quantitative ConsensusMap prior to this step.

	*/
	class OPENMS_DLLAPI ProteinInference
	{

		public:
	
			typedef Peak2D::IntensityType IntensityType;
			
			/// Constructor
			ProteinInference();
			
			/// copy constructor
	    ProteinInference(const ProteinInference& cp);
	
	    /// assignment operator
	    ProteinInference& operator = (const ProteinInference& rhs);
	
			/**
				@brief given a peptide quantitation, infer corresponding protein quantities
				
				Infers protein ratios from peptide ratios (currently using unique peptides only).
				Use the IDMapper class to add protein and peptide information to a 
				quantitative ConsensusMap prior to this step.
				
				@param consensus_map Peptide quantitation with ProteinIdentifications attached, where
							 Protein quantitation will be attached
				@param reference_map Index of (iTRAQ) reference channel within the consensus map
							 
				@throws Exception::MissingInformation if Protein/PeptideIdentifications are missing
			*/
			void infer(ConsensusMap& consensus_map,
								 const UInt reference_map);
			
			
		protected:
			
			void infer_(ConsensusMap& consensus_map, 
									const size_t protein_idenfication_index, 
									const UInt reference_map);
			
			bool sortByUnique_(std::vector< PeptideHit >& peptide_hits_local, const bool is_higher_score_better );
		
	}; // !class

} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_PROTEININFERENCE_H
