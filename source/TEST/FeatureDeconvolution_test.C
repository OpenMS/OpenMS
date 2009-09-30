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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

#include <coin/CoinMessageHandler.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureDeconvolution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureDeconvolution* ptr = 0;
START_SECTION(FeatureDeconvolution())
	ptr = new FeatureDeconvolution();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~FeatureDeconvolution())
	delete ptr;
END_SECTION

START_SECTION(void updateMembers_())
	NOT_TESTABLE
END_SECTION

START_SECTION(FeatureDeconvolution(const FeatureDeconvolution &source))
	FeatureDeconvolution fd;
	Param p;
	p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
	fd.setParameters(p);
	FeatureDeconvolution fd2(fd);
	FeatureDeconvolution fd_untouched;
	
	TEST_EQUAL(fd2.getParameters(), fd.getParameters())
	TEST_NOT_EQUAL(fd2.getParameters(), fd_untouched.getParameters())

END_SECTION

START_SECTION(FeatureDeconvolution& operator=(const FeatureDeconvolution &source))
	FeatureDeconvolution fd;
	Param p;
	p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
	fd.setParameters(p);
	FeatureDeconvolution fd2 = fd;
	FeatureDeconvolution fd_untouched;
	
	TEST_EQUAL(fd2.getParameters(), fd.getParameters())
	TEST_NOT_EQUAL(fd2.getParameters(), fd_untouched.getParameters())
END_SECTION


START_SECTION(void compute(const FeatureMapType &fm_in, FeatureMapType &fm_out, ConsensusMap &cons_map, ConsensusMap &cons_map_p))
//_CrtSetDbgFlag(_CrtSetDbgFlag(0)|_CRTDBG_CHECK_ALWAYS_DF);

	FeatureDeconvolution fd;
	FeatureMap<> fm_in, fm_out;
	ConsensusMap cm, cm2;
	FeatureXMLFile fl;
	fl.load(OPENMS_GET_TEST_DATA_PATH("FeatureDeconvolution_easy_input.featureXML"), fm_in);
	fd.compute(fm_in, fm_out, cm, cm2);

	String out_file;
	NEW_TMP_FILE(out_file)
	ConsensusXMLFile c1;

	c1.store(out_file,cm);

	FuzzyStringComparator fsc;
	fsc.setWhitelist (StringList::create("xml-stylesheet,map id"));
	fsc.setAcceptableAbsolute(0.01);
	bool cmp = fsc.compareFiles(out_file, OPENMS_GET_TEST_DATA_PATH("FeatureDeconvolution_easy_output.consensusXML"));
	TEST_EQUAL(cmp, true);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


