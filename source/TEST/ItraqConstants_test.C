// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
ItraqConstants* nullPointer = 0;
START_SECTION(ItraqConstants())
{
	ptr = new ItraqConstants();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
  ic.resize(3);
	ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
	ic[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
  ic[2].setMatrix<6,4>(ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX);

  {
		StringList ics = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::FOURPLEX, ic);
		StringList t_ics = StringList::create("114:0/1/5.9/0.2,115:0/2/5.6/0.1,116:0/3/4.5/0.1,117:0.1/4/3.5/0.1");
		TEST_EQUAL(ics, t_ics);
	}
	{	
		StringList ics = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::EIGHTPLEX, ic);
		StringList t_ics = StringList::create("113:0/0/6.89/0.22,114:0/0.94/5.9/0.16,115:0/1.88/4.9/0.1,116:0/2.82/3.9/0.07,117:0.06/3.77/2.99/0,118:0.09/4.71/1.88/0,119:0.14/5.66/0.87/0,121:0.27/7.44/0.18/0");
		TEST_EQUAL(ics, t_ics);
	}	 
  {
    StringList ics = ItraqConstants::getIsotopeMatrixAsStringList(ItraqConstants::TMT_SIXPLEX, ic);
    StringList t_ics = StringList::create("126:0/0/0/0,127:0/0/0/0,128:0/0/0/0,129:0/0/0/0,130:0/0/0/0,131:0/0/0/0");
    TEST_EQUAL(ics, t_ics);
  }
}
END_SECTION

START_SECTION((static void updateIsotopeMatrixFromStringList(const int itraq_type, const StringList &channels, IsotopeMatrices &isotope_corrections)))
{
	ItraqConstants::IsotopeMatrices ic;
  ic.resize(3);
	ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
	ic[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
  ic[2].setMatrix<6,4>(ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX);

//StringList t_ics = StringList::create("114:0/1/5.9/0.2,115:0/2/5.6/0.1,116:0/3/4.5/0.1,117:0.1/4/3.5/0.1"); // the default
	StringList t_ics = StringList::create("114:0/1/5.9/4.2,115:3/2/5.6/0.1,116:0/3/4.5/0.1,117:0.1/4/3.5/2");
	
	ic[0].setValue(0,3,4.2);
	ic[0].setValue(1,0,3); 
	ic[0].setValue(3,3,2);
	 
	ItraqConstants::IsotopeMatrices ic_new;
  ItraqConstants::updateIsotopeMatrixFromStringList(ItraqConstants::FOURPLEX, t_ics, ic_new);

  TEST_EQUAL(ic_new.size(), ic.size())
  for(ItraqConstants::IsotopeMatrices::size_type i = 0; i < ic_new.size() && i < ic.size(); ++i)
  {
    TEST_EQUAL(ic_new[i], ic[i])
  }

  // reset previously updated and update TMT isotope corrections
  ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
  ic[2].setValue(0,2,3.4);
  ic[2].setValue(1,0,2.1);
  ic[2].setValue(4,3,5.1);

  // StringList tmt_ics = StringList::create("126:0/0/0/0,127:0/0/0/0,128:0/0/0/0,129:0/0/0/0,130:0/0/0/0,131:0/0/0/0"); // the original one
  StringList tmt_ics = StringList::create("126:0/0/3.4/0,127:2.1/0/0/0,128:0/0/0/0,129:0/0/0/0,130:0/0/0/5.1,131:0/0/0/0");

  ItraqConstants::IsotopeMatrices ic_tmt;
  ItraqConstants::updateIsotopeMatrixFromStringList(ItraqConstants::TMT_SIXPLEX, tmt_ics, ic_tmt);

  TEST_EQUAL(ic_new.size(), ic.size())
  for(ItraqConstants::IsotopeMatrices::size_type i = 0; i < ic_tmt.size() && i < ic.size(); ++i)
  {
    TEST_EQUAL(ic_tmt[i], ic[i])
  }
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
  
  ItraqConstants::ChannelMapType mapTMT;
  ItraqConstants::initChannelMap(ItraqConstants::TMT_SIXPLEX, mapTMT);

  TEST_EQUAL(6, mapTMT.size());
  TEST_EQUAL(mapTMT[126].id, 0);
  TEST_EQUAL(mapTMT[126].active, false);
  TEST_EQUAL(mapTMT[129].id, 3);
  TEST_EQUAL(mapTMT[129].active, false);
}
END_SECTION

START_SECTION((static void updateChannelMap(const StringList& active_channels, ChannelMapType& map)))
{
	StringList active_channels = StringList::create("114:myReference");
	ItraqConstants::ChannelMapType map;
	
	ItraqConstants::initChannelMap(ItraqConstants::FOURPLEX, map);
	
  ItraqConstants::updateChannelMap(active_channels, map);
  
  TEST_EQUAL(map[114].description, String("myReference"))
	TEST_EQUAL(map[114].active, true);
  
  // TMT
  StringList activeTmtChannels = StringList::create("126:myReference,129:treated,131:control");
  ItraqConstants::ChannelMapType tmtMap;

  ItraqConstants::initChannelMap(ItraqConstants::TMT_SIXPLEX, tmtMap);

  ItraqConstants::updateChannelMap(activeTmtChannels, tmtMap);

  TEST_EQUAL(tmtMap[126].description, String("myReference"))
  TEST_EQUAL(tmtMap[126].active, true);
  TEST_EQUAL(tmtMap[127].description, String(""))
  TEST_EQUAL(tmtMap[127].active, false);
  TEST_EQUAL(tmtMap[128].description, String(""))
  TEST_EQUAL(tmtMap[128].active, false);
  TEST_EQUAL(tmtMap[129].description, String("treated"))
  TEST_EQUAL(tmtMap[129].active, true);
  TEST_EQUAL(tmtMap[130].description, String(""))
  TEST_EQUAL(tmtMap[130].active, false);
  TEST_EQUAL(tmtMap[131].description, String("control"))
  TEST_EQUAL(tmtMap[131].active, true);
}
END_SECTION

START_SECTION((static Matrix<double> translateIsotopeMatrix(const int &itraq_type, const IsotopeMatrices &isotope_corrections)))
{
  ItraqConstants::IsotopeMatrices ic;
  ic.resize(3);
	ic[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
	ic[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
  ic[2].setMatrix<6,4>(ItraqConstants::ISOTOPECORRECTIONS_TMT_SIXPLEX);

  Matrix<double> channel_frequency = ItraqConstants::translateIsotopeMatrix(ItraqConstants::FOURPLEX, ic);
	
	std::cout << "CF: \n" << channel_frequency << "\n";
	TEST_REAL_SIMILAR(channel_frequency.getValue(0,0), 0.929)
	TEST_REAL_SIMILAR(channel_frequency.getValue(3,0), 0)

  channel_frequency = ItraqConstants::translateIsotopeMatrix(ItraqConstants::EIGHTPLEX, ic);

  std::cout << "CF: \n" << channel_frequency << "\n";
  /*
    0.9289 0.0094      0      0      0      0      0      0
    0.0689   0.93 0.0188      0      0      0      0      0
    0.0022  0.059 0.9312 0.0282 0.0006      0      0      0
         0 0.0016  0.049 0.9321 0.0377 0.0009      0      0
         0      0  0.001  0.039 0.9329 0.0471 0.0014      0
         0      0      0 0.0007 0.0288 0.9332 0.0566      0
         0      0      0      0      0 0.0188 0.9333 0.0027
         0      0      0      0      0      0      0 0.9211

  */
  // test lower right triangle
  TEST_REAL_SIMILAR(channel_frequency.getValue(6,7), 0.0027)
  TEST_REAL_SIMILAR(channel_frequency.getValue(7,7), 0.9211)
  TEST_REAL_SIMILAR(channel_frequency.getValue(7,6), 0.0000)

  channel_frequency = ItraqConstants::translateIsotopeMatrix(ItraqConstants::TMT_SIXPLEX, ic);
  std::cout << "CF: \n" << channel_frequency << "\n";
  TEST_REAL_SIMILAR(channel_frequency.getValue(0,0), 1.0)
  TEST_REAL_SIMILAR(channel_frequency.getValue(1,0), 0.0)
  TEST_REAL_SIMILAR(channel_frequency.getValue(0,1), 0.0)
  TEST_REAL_SIMILAR(channel_frequency.getValue(3,3), 1.0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



