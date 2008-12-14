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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ITRAQQUANTIFIER_H
#define OPENMS_ANALYSIS_QUANTITATION_ITRAQQUANTIFIER_H

#include <vector>

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>


namespace OpenMS
{

	/**
		@brief
		
		@todo Docu (Chris)
	*/
	class OPENMS_DLLAPI ItraqQuantifier
		: public DefaultParamHandler,
			public ItraqConstants
	{

	public:

		typedef ItraqConstants::ChannelInfo ChannelInfo;
		typedef ItraqConstants::ChannelMapType ChannelMapType;

		
		/// default isotope correction factors
		static const double ISOTOPECORRECTIONS_FOURPLEX[4][4];
		static const double ISOTOPECORRECTIONS_EIGHTPLEX[8][4];
		static const Int CHANNEL_COUNT[];
		
		/// Constructor (assuming 4-plex experiment)
		ItraqQuantifier();
		
		/// Constructor with iTRAQ-type (either ItraqQuantifier::FOURPLEX or ItraqQuantifier::EIGHTPLEX)
		ItraqQuantifier(Int itraq_type);

		/// Constructor with iTRAQ-type (either ItraqQuantifier::FOURPLEX or ItraqQuantifier::EIGHTPLEX) and Param
		ItraqQuantifier(Int itraq_type, const Param& param);

		/// copy constructor
    ItraqQuantifier(const ItraqQuantifier& cp);

    /// assignment operator
    ItraqQuantifier& operator = (const ItraqQuantifier& rhs);

		/**
		 *	@brief using the raw iTRAQ intensities we apply isotope correction, normalization (using median) and protein inference
		 *	
		 *	@param consensus_map_in Raw iTRAQ intensities from previous step
		 *	@param consensus_map_out Postprocessed iTRAQ ratios for peptides
		 *
		 *	@throws Exception::FailedAPICall is least-squares fit fails
		 *	@throws Exception::InvalidParameter if parameter is invalid (e.g. reference_channel)
		 */
		void run(const ConsensusMap& consensus_map_in, 
						 ConsensusMap& consensus_map_out
						 );
		
		/**
		 *	@brief using the raw iTRAQ intensities we apply isotope correction, normalization (using median) and protein inference
		 *	
		 *	@param consensus_map_in Raw iTRAQ intensities from previous step
		 *	@param peptide_ids List of peptides identified by a search engine on the same MS² dataset
		 *	@param protein_ids List of proteins inferred from peptides
		 *	@param consensus_map_out Postprocessed iTRAQ ratios for Proteins (if provided) or Peptides otherwise
		 *
		 *	@throws Exception::FailedAPICall is least-squares fit fails
		 *	@throws Exception::InvalidParameter if parameter is invalid (e.g. reference_channel)
		 */
		void run(const ConsensusMap& consensus_map_in, 
						 const std::vector< PeptideIdentification > &peptide_ids,
						 const std::vector< ProteinIdentification > &protein_ids, 
						 ConsensusMap& consensus_map_out
						 );

	protected:
		
		void setDefaultParams_();
		
		void updateMembers_();
		
	private:
		
		/// initialize
		void initIsotopeCorrections_();
		
		void reconstruct_channel_info_(const ConsensusMap& consensus_map);
			
		/// either ItraqQuantifier::FOURPLEX or ItraqQuantifier::EIGHTPLEX
		Int itraq_type_;
		
		/// map the channel-name (eg 114) onto its channel_info
		/// the channel-description is also the id-string in the mapList section of the ConsensusMap
		ChannelMapType channel_map_;	

		/// Matrixes with isotope correction values (one for each plex-type)
		std::vector< Matrix<double> > isotope_corrections_;
		
	}; // !class

} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ITRAQQUANTIFIER_H
