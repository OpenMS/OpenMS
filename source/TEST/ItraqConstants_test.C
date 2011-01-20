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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ItraqConstants, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqConstants* ptr = 0;
START_SECTION(ItraqConstants())
{
	ptr = new ItraqConstants();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ItraqConstants())
{
	delete ptr;
}
END_SECTION

START_SECTION((static StringList getIsotopeMatrixAsStringList(const int itraq_type, const IsotopeMatrices &isotope_corrections)))
{
	ItraqConstants::IsotopeMatrices ic;
	ic.resize(2);
	ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
	ic[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);

	{
		StringList ics = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::FOURPLEX, ic);
		StringList t_ics = StringList::create("114:0/1/5.9/0.2,115:0/2/5.6/0.1,116:0/3/4.5/0.1,117:0.1/4/3.5/0");
		TEST_EQUAL(ics, t_ics);
	}
	{	
		StringList ics = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::EIGHTPLEX, ic);
		StringList t_ics = StringList::create("113:0/2.5/6.89/0.22,114:0/0.94/5.9/0.16,115:0/1.88/4.9/0.1,116:0/2.82/3.9/0.07,117:0.06/3.77/2.88/0,118:0.09/4.71/1.88/0,119:0.14/5.66/0.87/0,121:0.27/7.44/0.18/0");
		TEST_EQUAL(ics, t_ics);
	}	 
}
END_SECTION

START_SECTION((static void updateIsotopeMatrixFromStringList(const int itraq_type, const StringList &channels, IsotopeMatrices &isotope_corrections)))
{
	ItraqConstants::IsotopeMatrices ic;
	ic.resize(2);
	ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
	ic[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);

//StringList t_ics = StringList::create("114:0/1/5.9/0.2,115:0/2/5.6/0.1,116:0/3/4.5/0.1,117:0.1/4/3.5/0"); // the default
	StringList t_ics = StringList::create("114:0/1/5.9/4.2,115:3/2/5.6/0.1,116:0/3/4.5/0.1,117:0.1/4/3.5/2");
	
	ic[0].setValue(0,3,4.2);
	ic[0].setValue(1,0,3); 
	ic[0].setValue(3,3,2);
	 
	ItraqConstants::IsotopeMatrices ic_new;
  ItraqConstants::updateIsotopeMatrixFromStringList(ItraqConstants::FOURPLEX, t_ics, ic_new);
  
  TEST_EQUAL(ic_new[0], ic[0])
  
}
END_SECTION

START_SECTION((static void initChannelMap(const int itraq_type, ChannelMapType &map)))
{
	ItraqConstants::ChannelMapType map;
  ItraqConstants::initChannelMap(ItraqConstants::EIGHTPLEX, map);
  
  TEST_EQUAL(8, map.size());
  TEST_EQUAL(map[119].id, 6);
  TEST_EQUAL(map[119].active, false);
  
	ItraqConstants::ChannelMapType map4;
  ItraqConstants::initChannelMap(ItraqConstants::FOURPLEX, map4);
  
  TEST_EQUAL(4, map4.size());
  TEST_EQUAL(map4[114].id, 0);
  TEST_EQUAL(map4[114].active, false);
  
  
}
END_SECTION

START_SECTION((static void updateChannelMap(StringList active_channels, ChannelMapType &map)))
{
	StringList active_channels = StringList::create("114:myReference");
	ItraqConstants::ChannelMapType map;
	
	ItraqConstants::initChannelMap(ItraqConstants::FOURPLEX, map);
	
  ItraqConstants::updateChannelMap(active_channels, map);
  
  TEST_EQUAL(map[114].description, String("myReference"))
	TEST_EQUAL(map[114].active, true);
  
}
END_SECTION

START_SECTION((static Matrix<double> translateIsotopeMatrix(const int &itraq_type, const IsotopeMatrices &isotope_corrections)))
{
  ItraqConstants::IsotopeMatrices ic;
	ic.resize(2);
	ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
	ic[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
  Matrix<double> channel_frequency = ItraqConstants::translateIsotopeMatrix(ItraqConstants::FOURPLEX, ic);
	
	std::cout << "CF: " << channel_frequency << "\n";
	TEST_REAL_SIMILAR(channel_frequency.getValue(0,0), 0.929)
	TEST_REAL_SIMILAR(channel_frequency.getValue(3,0), 0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



