// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
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
			Peak2D::CoordinateType center; // expected centroid of peak in MZ
			bool active; // channel actually added to the experiment?
		};

		/// maps iTRAQ channel (e.g. 117) to more information
		typedef Map<Int, ChannelInfo > ChannelMapType;

		/// (user defined?) isotope correction matrix in (-2, -1, +1, +2) row style
		typedef std::vector< Matrix<double> > IsotopeMatrices;
		
		/// channel names for 4plex( 114, 115, 116, 117)
		static const Int CHANNELS_FOURPLEX[4][1];
		/// channel names for 8plex( 113, 114, 115, 116, 117, 118, 119, 121)
		static const Int CHANNELS_EIGHTPLEX[8][1];

		/// default isotope correction matrix (4 plex)
		static const double ISOTOPECORRECTIONS_FOURPLEX[4][4];
		/// default isotope correction matrix (8 plex)
		static const double ISOTOPECORRECTIONS_EIGHTPLEX[8][4];

		/**
			@brief convert isotope correction matrix to stringlist

			Each line is converted into a string of the format &lt;channel&gt;:&lt;-2Da&gt;/&lt;-1Da&gt;/&lt;+1Da&gt;/&lt;+2Da&gt; ; e.g. '114:0/0.3/4/0' 
			Useful for creating parameters or debug output.

			@param itraq_type Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
			@param isotope_corrections Vector of the two matrices (4plex, 8plex).
		**/
		static StringList getIsotopeMatrixAsStringList(const int itraq_type, const IsotopeMatrices& isotope_corrections);

		/**
			@brief convert strings to isotope correction matrix rows

			Each string of format &lt;channel&gt;:&lt;-2Da&gt;/&lt;-1Da&gt;/&lt;+1Da&gt;/&lt;+2Da&gt; ; e.g. '114:0/0.3/4/0' 
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
			Each input string must have the format &lt;channel&gt;:&lt;description&gt;, e.g. "114:myref","115:liver"

			@param active_channels StringList with channel and description
			@param map Storage to update
		**/
		static void updateChannelMap(const StringList& active_channels, ChannelMapType& map);

		/**
			@brief translate isotope correction matrix in -2,-1,+1,+2 form into 114,115,116,117 format
			
			Translates e.g. ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX matrix into a 8x8 matrix which
			maps how channel (row) distributes its tags onto other channels (columns).
			
			@param itraq_type Should be of values from enum ITRAQ_TYPES
			@param isotope_corrections isotope correction matrix in -2...+2 form
		**/
		static Matrix<double> translateIsotopeMatrix(const int& itraq_type, const IsotopeMatrices& isotope_corrections);
		
	}; // !class
	
} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ITRAQCONSTANTS_H
 
