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
// $Maintainer: Alexandra Scherbart $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/PIP/LocalLinearMap.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(LocalLinearMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LocalLinearMap* ptr = 0;
LocalLinearMap llm;


START_SECTION(LocalLinearMap())
	ptr = new LocalLinearMap();
	TEST_NOT_EQUAL(ptr, 0)
	TEST_EQUAL(ptr->getLLMParam().xdim, llm.getLLMParam().xdim)
	TEST_EQUAL(ptr->getLLMParam().ydim, llm.getLLMParam().ydim)
	TEST_EQUAL(ptr->getLLMParam().radius, llm.getLLMParam().radius)
END_SECTION

START_SECTION((virtual ~LocalLinearMap()))
	delete ptr;
END_SECTION

START_SECTION(const LLMParam& getLLMParam() const)
	TEST_EQUAL(llm.getLLMParam().xdim, 1)
	TEST_EQUAL(llm.getLLMParam().ydim, 2)
	TEST_EQUAL(llm.getLLMParam().radius, 0.4)
END_SECTION

START_SECTION(const Matrix<DoubleReal>& getCodebooks() const)
	TEST_EQUAL(llm.getCodebooks().rows(), 2)
	TEST_EQUAL(llm.getCodebooks().cols(), 18)
	//0.0163859	0.9420950
	//0.0368383	-0.4910166
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(0,0), 0.030113)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(0,1), 0.01550)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(1,0), 0.0)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(1,1), 0.0)
END_SECTION

START_SECTION(const Matrix<DoubleReal>& getMatrixA() const)
	TEST_EQUAL(llm.getMatrixA().rows(), 2)
	TEST_EQUAL(llm.getMatrixA().cols(), 18)
	//-0.4431946	0.2819091
	//-0.5988132	-0.1837768
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(0,0), 3.31028978)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(0,1), 0.0)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(1,0), 0.0)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(1,1), 0.0)
END_SECTION

START_SECTION(const vector<DoubleReal>& getVectorWout() const)
	TEST_EQUAL(llm.getVectorWout().size(), 2)
	//4.205033
	//4.205731
	TEST_REAL_SIMILAR(llm.getVectorWout()[0], 3.8171745)
	TEST_REAL_SIMILAR(llm.getVectorWout()[1], 0.0)
END_SECTION

START_SECTION(const Matrix<DoubleReal>& getCord() const)
	TEST_EQUAL(llm.getCord().rows(), 2)
	TEST_EQUAL(llm.getCord().cols(), 2)
	TEST_REAL_SIMILAR(llm.getCord().getValue(0,0), 0)
	TEST_REAL_SIMILAR(llm.getCord().getValue(0,1), 0)
	TEST_REAL_SIMILAR(llm.getCord().getValue(1,0), 0)
	TEST_REAL_SIMILAR(llm.getCord().getValue(1,1), 1)
END_SECTION

START_SECTION(vector<DoubleReal> neigh(const Matrix<UInt>& rhs1, UInt rhs2, DoubleReal rhs3))
	vector<DoubleReal> nei1 = llm.neigh(llm.getCord(), 0, llm.getLLMParam().radius);
	TEST_EQUAL(nei1[0], 1)
	TEST_REAL_SIMILAR(nei1[1], 0.04393693) 	
END_SECTION

START_SECTION(void normalizeVector(std::vector<DoubleReal>& rhs))
	NOT_TESTABLE
END_SECTION


END_TEST
