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
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(MapAlignmentEvaluationAlgorithmPrecision, "$Id MapAlignmentEvaluationAlgorithmPrecision_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignmentEvaluationAlgorithmPrecision* ptr = 0;
START_SECTION((MapAlignmentEvaluationAlgorithmPrecision()))
	ptr = new MapAlignmentEvaluationAlgorithmPrecision();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~MapAlignmentEvaluationAlgorithmPrecision()))
	delete ptr;
END_SECTION

START_SECTION((static MapAlignmentEvaluationAlgorithm* create()))
	MapAlignmentEvaluationAlgorithm* ptr2 = 0;
	ptr2 = MapAlignmentEvaluationAlgorithmPrecision::create();
	TEST_NOT_EQUAL(ptr2, 0)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(MapAlignmentEvaluationAlgorithmPrecision::getProductName(),"precision")
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap &consensus_map_in, const ConsensusMap &consensus_map_gt, const DoubleReal &rt_dev, const DoubleReal &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge, DoubleReal &out)))
	MapAlignmentEvaluationAlgorithmPrecision maea;
	ConsensusMap in;
	ConsensusMap gt;
	DoubleReal out; 
	
	ConsensusXMLFile consensus_xml_file_in;
	consensus_xml_file_in.load( OPENMS_GET_TEST_DATA_PATH("MapAlignmentEvaluationAlgorithm_in.consensusXML"), in );
		
	ConsensusXMLFile consensus_xml_file_gt;
	consensus_xml_file_gt.load( OPENMS_GET_TEST_DATA_PATH("MapAlignmentEvaluationAlgorithm_gt.consensusXML"), gt );
	
	maea.evaluate(in, gt, 0.1, 0.1, 100, true, out);

	TEST_REAL_SIMILAR(out, 0.757143)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

