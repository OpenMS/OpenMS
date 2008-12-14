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

#ifndef OPENMS_ANALYSIS_QUANTITATION_CHANNELEXTRACTOR_H
#define OPENMS_ANALYSIS_QUANTITATION_CHANNELEXTRACTOR_H

#include <vector>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>

namespace OpenMS
{
	
	/**
		@brief
		
		@todo Docu (Chris)
	*/
	class OPENMS_DLLAPI ItraqChannelExtractor
		: public DefaultParamHandler,
			public ItraqConstants
	{

	public:
		
		typedef ItraqConstants::ChannelMapType ChannelMapType;

		/// Constructor (assuming 4plex)
		ItraqChannelExtractor();

		/// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES)
		ItraqChannelExtractor(Int itraq_type);
		
		/// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES) and param
		ItraqChannelExtractor(Int itraq_type, const Param& param);

		/// copy constructor
		ItraqChannelExtractor(const ItraqChannelExtractor& cp);

		/// assignment operator
		ItraqChannelExtractor& operator = (const ItraqChannelExtractor& rhs);

		
		/// @brief extracts the iTRAQ channels from the MS data and stores intensity values in a consensus map
		///
		/// @param ms_exp_data Raw data to read
		/// @param consensus_map 
		void run(const MSExperiment<Peak1D>& ms_exp_data, ConsensusMap& consensus_map);

	protected:
		
		void setDefaultParams_();
		
		/// implemented for DefaultParamHandler
		void updateMembers_();
			
		
	private:
		
		/// initialize
		void init_();
		
		/// set to either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX
		Int itraq_type_;
		
		/// map the channel-name (eg 114) onto its description and the centroid mass
		/// the channel-name is also the id-string in the mapList section of the ConsensusMap
		ChannelMapType channel_map_;	
		
		/// valid names for 4 and 8plex, ie 114,115,116,117 for 4plex
		std::vector< Matrix<Int> > channel_names_;

	}; // !class

	
} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_CHANNELEXTRACTOR_H
 
