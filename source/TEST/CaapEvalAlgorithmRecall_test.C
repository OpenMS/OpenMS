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
// $Maintainer: Katharina Albers, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithmRecall.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(CaapEvalAlgorithmRecall, "$Id CaapEvalAlgorithmRecall_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CaapEvalAlgorithmRecall* ptr = 0;
START_SECTION((CaapEvalAlgorithmRecall()))
	ptr = new CaapEvalAlgorithmRecall();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~CaapEvalAlgorithmRecall()))
	delete ptr;
END_SECTION

START_SECTION((static CaapEvalAlgorithmRecall* create()))
	CaapEvalAlgorithm* ptr2 = 0;
	ptr2 = CaapEvalAlgorithmRecall::create();
	TEST_NOT_EQUAL(ptr2, 0)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(CaapEvalAlgorithmRecall::getProductName(),"recall")
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap& mapin1, const ConsensusMap& mapin2, DoubleReal& out)))
	CaapEvalAlgorithmRecall cae;
	ConsensusMap in;
	ConsensusMap gt;
	DoubleReal out; 
	
	ConsensusXMLFile consensus_xml_file_in;
	consensus_xml_file_in.load( OPENMS_GET_TEST_DATA_PATH("CaapEvalAlgorithm_in.consensusXML"), in );
		
	ConsensusXMLFile consensus_xml_file_gt;
	consensus_xml_file_gt.load( OPENMS_GET_TEST_DATA_PATH("CaapEvalAlgorithm_gt.consensusXML"), gt );
	
	cae.evaluate(in, gt, out);

	TEST_REAL_SIMILAR(out, 0.5)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
