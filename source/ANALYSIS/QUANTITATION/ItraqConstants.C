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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>


namespace OpenMS {

	// number of channels for iTRAQ types. (make sure it corresponds to enum ITRAQ_TYPES)
	const Int ItraqConstants::CHANNEL_COUNT[2] = {4, 8};

	const Int ItraqConstants::CHANNELS_FOURPLEX[4][1] = {{114}, {115}, {116}, {117}};
	const Int ItraqConstants::CHANNELS_EIGHTPLEX[8][1] = {{113}, {114}, {115}, {116}, {117}, {118}, {119}, {121}};

	// currently from http://www.matrixscience.com/help/quant_config_help.html
	const double ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX[4][4] = {
		{0.0, 1.0, 5.9, 0.2},		//114
		{0.0, 2.0, 5.6, 0.1},
		{0.0, 3.0, 4.5, 0.1},
		{0.1, 4.0, 3.5, 0.0}		//117
	};
	
	//taken from Applied Biosystems Website
	// http://faqs.appliedbiosystems.com/cgi-bin/appliedbio.cfg/php/enduser/std_adp.php?p_faqid=3671
	const double ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX[8][4] = {
		{0.00, 2.50, 6.89, 0.22},		//113
		{0.00, 0.94, 5.90, 0.16},
		{0.00, 1.88, 4.90, 0.10},
		{0.00, 2.82, 3.90, 0.07},
		{0.06, 3.77, 2.88, 0.00},
		{0.09, 4.71, 1.88, 0.00},
		{0.14, 5.66, 0.87, 0.00},
		{0.27, 7.44, 0.18, 0.00}		//121
	};

	StringList ItraqConstants::getIsotopeMatrixAsStringList(const int itraq_type, const IsotopeMatrices& isotope_corrections)
	{
		OPENMS_PRECONDITION(itraq_type < SIZE_OF_ITRAQ_TYPES && itraq_type>=0, "Error while trying to access invalid isotope correction matrix.");
		
		StringList isotopes;
		std::vector< Matrix<Int> > channel_names(2);
		channel_names[0].setMatrix<4,1>(CHANNELS_FOURPLEX);
		channel_names[1].setMatrix<8,1>(CHANNELS_EIGHTPLEX);
		for (Int i=0; i<CHANNEL_COUNT[itraq_type]; ++i)
		{
			String line = String(channel_names[itraq_type].getValue(i,0)) + ":";
			for (Int j=0;j<3;++j)
			{
				line += String(isotope_corrections[itraq_type].getValue(i,j)) + "/";
			}
			line += String(isotope_corrections[itraq_type].getValue(i,3));
			isotopes.push_back(line);
		} 

		return isotopes;
	}

	void ItraqConstants::updateIsotopeMatrixFromStringList(const int itraq_type, const StringList& channels, IsotopeMatrices& isotope_corrections)
	{	
		isotope_corrections.resize(2);
		isotope_corrections[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
		isotope_corrections[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
		
		// split the channels key:name pairs apart
		for (StringList::const_iterator it=channels.begin();it!=channels.end();++it)
		{
			StringList result;
			it->split(':',result);
			if (result.size()!=2)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; expected one ':', got this: '" + (*it) + "'");
			}
			result[0] = result[0].trim(); // hold channel name
			result[1] = result[1].trim(); // holds 4 values

			Int channel = result[0].toInt();
			Int line = (itraq_type == FOURPLEX ? channel-114 : channel-113);
			if ((itraq_type == FOURPLEX && (line<0 || line>3))
					|| 
					((itraq_type == EIGHTPLEX && (line<0 || line>8)) || channel==120))
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,String("ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; channel-name is not valid for ") + String(itraq_type==FOURPLEX ? "4plex": "8plex") + String(": '") + result[0] + String("'"));
			}
			// if set to 121 we still want to change line 7 of the matrix
			if (line==8) line=7;

			StringList corrections;
			result[1].split('/',corrections);
			if (corrections.size()!=4)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqQuantifier: Invalid entry in Param 'isotope_correction_values'; expected four correction values separated by '&', got this: '" + result[1] + "'");
			}

			// overwrite line in Matrix with custom values
			isotope_corrections[itraq_type].setValue(line,0, corrections[0].toDouble());
			isotope_corrections[itraq_type].setValue(line,1, corrections[1].toDouble());
			isotope_corrections[itraq_type].setValue(line,2, corrections[2].toDouble());
			isotope_corrections[itraq_type].setValue(line,3, corrections[3].toDouble());

			#ifdef ITRAQ_DEBUG
			std::cout << "Channel " << channel << " has values " << corrections << std::endl;
			#endif
		}
	}

	void ItraqConstants::initChannelMap(const int itraq_type, ChannelMapType& map)
	{
		/// valid names for 4 and 8plex, ie 114,115,116,117 for 4plex
		std::vector< Matrix<Int> > channel_names;
		channel_names.resize(2);
		channel_names[0].setMatrix<4,1>(CHANNELS_FOURPLEX);
		channel_names[1].setMatrix<8,1>(CHANNELS_EIGHTPLEX);

		map.clear();
		for (Size i=0; i < channel_names[itraq_type].rows(); ++i)
		{
			ChannelInfo info;
			info.description = "";
			info.name = channel_names[itraq_type].getValue(i,0);
			info.id = (Int)i;
			info.center = double(info.name) + 0.11;
			info.active = false;
			map[info.name] = info;
		}

		#ifdef ITRAQ_DEBUG
		std::cout << "INIT: map has " << map.size() << " entries!" << std::endl;
		#endif		
	}

	void ItraqConstants::updateChannelMap(StringList active_channels, ChannelMapType& map)
	{
		// split the channels key:name pairs apart
		for (StringList::const_iterator it=active_channels.begin();it!=active_channels.end();++it)
		{
			StringList result;
			it->split(':',result);
			if (result.size()!=2)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqConstants: Invalid entry in Param 'channel_active'; expected one semicolon ('" + (*it) + "')");
			}
			result[0] = result[0].trim();
			result[1] = result[1].trim();
			if (result[0]==String::EMPTY || result[1]==String::EMPTY)
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqConstants: Invalid entry in Param 'channel_active'; key or value is empty ('" + (*it) + "')");
			}
			Int channel = result[0].toInt();
			if (map.find(channel) == map.end())
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"ItraqConstants: Invalid entry in Param 'channel_active'; channel is not valid ('" + String(channel) + "')");
			}
			// update name (description) of channel
			map[channel].description = result[1];
			map[channel].active = true;

			#ifdef ITRAQ_DEBUG
			std::cout << "Channel " << channel << " has description " << channel_map_[channel].description << " and center " << channel_map_[channel].center << std::endl;
			#endif
		}
	}

	Matrix<double> ItraqConstants::translateIsotopeMatrix(const int& itraq_type, const IsotopeMatrices& isotope_corrections)
	{
		// translate isotope_corrections to a channel_frequency matrix
		Matrix<double> channel_frequency(CHANNEL_COUNT[itraq_type], CHANNEL_COUNT[itraq_type]);
		for (Int i=0; i < CHANNEL_COUNT[itraq_type]; ++i)
		{
			for (Int j=0; j < CHANNEL_COUNT[itraq_type]; ++j)
			{
				// diagonal (should be close to 1 = 100%)
				if (i==j)
				{
					double val = 1.0;

					// subtract all isotope deviations of row i
					for (Int col_idx=0; col_idx < 4; ++col_idx) 
					{
						val += -isotope_corrections[itraq_type].getValue(i,col_idx) / 100;
					}
					channel_frequency.setValue(i,j,val);
				}
				else
				{ // from mass i to mass j (directly copy the deviation)
					if (i-j<=2 && i-j>0)
					{
						channel_frequency.setValue(j,i, isotope_corrections[itraq_type].getValue(i,j-i+2) / 100);
					}
					else if (j-i<=2 && j-i>0)
					{
						channel_frequency.setValue(j,i, isotope_corrections[itraq_type].getValue(i,j-i+1) / 100);
					}
				}
			}
		}
		
		return channel_frequency;
	}

} // !namespace
