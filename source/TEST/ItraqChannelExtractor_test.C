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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ItraqChannelExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqChannelExtractor* ptr = 0;
START_SECTION(ItraqChannelExtractor())
{
	ptr = new ItraqChannelExtractor();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ItraqChannelExtractor())
{
	delete ptr;
}
END_SECTION

START_SECTION((ItraqChannelExtractor(Int itraq_type)))
{
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX);
	TEST_EQUAL((StringList) ice.getParameters().getValue("channel_active"), StringList::create("113:liver,117:lung"));
  ItraqChannelExtractor ice2(ItraqChannelExtractor::FOURPLEX);
	TEST_EQUAL((StringList) ice2.getParameters().getValue("channel_active"), StringList::create("114:liver,117:lung"));
}
END_SECTION

START_SECTION((ItraqChannelExtractor(Int itraq_type, const Param &param)))
{
	Param p;
	p.setValue("reporter_mass_shift", 0.1234);
	p.setValue("channel_active", StringList::create("121:this is a test"));
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX, p);
	TEST_EQUAL((double) ice.getParameters().getValue("reporter_mass_shift"), 0.1234);
	TEST_EQUAL((StringList) ice.getParameters().getValue("channel_active"), StringList::create("121:this is a test"));
	
	// this should go wrong
	p.setValue("channel_active", StringList::create("120:channel non existant"));	
	TEST_EXCEPTION(Exception::InvalidParameter, ItraqChannelExtractor ice2(ItraqChannelExtractor::EIGHTPLEX, p));	
}
END_SECTION

START_SECTION((ItraqChannelExtractor(const ItraqChannelExtractor &cp)))
{
	Param p;
	p.setValue("reporter_mass_shift", 0.1234);
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX, p);
	ItraqChannelExtractor ice_cp(ice);
	
	TEST_EQUAL(ice_cp.getParameters(), ice.getParameters());
}
END_SECTION

START_SECTION((ItraqChannelExtractor& operator=(const ItraqChannelExtractor &rhs)))
{
	Param p;
	p.setValue("reporter_mass_shift", 0.1234);
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX, p);
	ItraqChannelExtractor ice_cp;
	ice_cp=ice;
	
	TEST_EQUAL(ice_cp.getParameters(), ice.getParameters());
}
END_SECTION



START_SECTION((void run(const MSExperiment< Peak1D > &ms_exp_data, ConsensusMap &consensus_map)))
{
	MzDataFile mz_data_file;
	MSExperiment<Peak1D > exp;
	mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.mzData"),exp);
	Param p;
	p.setValue("channel_active", StringList::create("114:ref,115:something,116:else"));
	p.setValue("select_activation","");
  ItraqChannelExtractor ice(ItraqChannelExtractor::FOURPLEX, p);
	ConsensusMap cm_out;
	ice.run(exp, cm_out);
	
	ConsensusXMLFile cm_file;
	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);
	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.consensusXML"));
	
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



