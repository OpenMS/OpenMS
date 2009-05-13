// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_ANALYSIS_QUANTITATION_ITRAQCONSTANTS_H
#define OPENMS_ANALYSIS_QUANTITATION_ITRAQCONSTANTS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{

	/**
		@brief Some constants used throughout iTRAQ classes.
		
		Constants for iTRAQ experiments and a ChannelInfo structure to store information about a single channel.

	*/
	class OPENMS_DLLAPI ItraqConstants
	{

	public:

		static const Int CHANNEL_COUNT[];
		enum ITRAQ_TYPES {FOURPLEX=0, EIGHTPLEX, SIZE_OF_ITRAQ_TYPES};
		
		/// stores information on an iTRAQ channel
		struct ChannelInfo
		{
			String description; // description given by experimentator (e.g. lung tissue)
			Int name; // 114-117 or 113 to 121
			Int id;		// 0-4 or 0-8
			Peak2D::CoordinateType center; // expected centoid of peak in MZ
			bool active; // channel actually added to the experiment?
		};

		typedef Map<Int, ChannelInfo > ChannelMapType;
		typedef std::vector< Matrix<double> > IsotopeMatrices;
		
		static const Int CHANNELS_FOURPLEX[4][1];
		static const Int CHANNELS_EIGHTPLEX[8][1];

		/// default isotope correction matrices
		static const double ISOTOPECORRECTIONS_FOURPLEX[4][4];
		static const double ISOTOPECORRECTIONS_EIGHTPLEX[8][4];

		/**
			@brief convert isotope correction matrix to stringlist

			Each line is converted into a string of the format <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' 
			Useful for creating parameters or debug output.

			@param itraq_type Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
			@param isotope_corrections Vector of the two matrices (4plex, 8plex).
		**/
		static StringList getIsotopeMatrixAsStringList(const int itraq_type, const IsotopeMatrices& isotope_corrections);

		/**
			@brief convert strings to isotope correction matrix rows

			Each string of format <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' 
			is parsed and the corresponding channel(row) in the matrix is updated.
			Not all channels need to be present, missing channels will be left untouched.
			Useful to update the matrix with user isotope correction values.

			@param itraq_type Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
			@param channels New channel isotope values as strings
			@param isotope_corrections Vector of the two matrices (4plex, 8plex).
		**/
		static void updateIsotopeMatrixFromStringList(const int itraq_type, const StringList& channels, IsotopeMatrices& isotope_corrections);

		/**
			@brief information about an iTRAQ channel

			State, name and expected mz-position of iTRAQ channels are initialized.

			@param itraq_type Should be of values from enum ITRAQ_TYPES
			@param map Storage to initialize
		**/
		static void initChannelMap(const int itraq_type, ChannelMapType& map);

		/**
			@brief activate & annotate channels

			State and description of iTRAQ channels are updated.
			Each input string must have the format <channel>:<description>, e.g. "114:myref","115:liver"

			@param active_channels StringList with channel and description
			@param map Storage to update
		**/
		static void updateChannelMap(StringList active_channels, ChannelMapType& map);

	}; // !class
	
} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ITRAQCONSTANTS_H
 
