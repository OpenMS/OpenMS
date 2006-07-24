// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:Cornelia Friedle $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/KERNEL/DSpectrum.h>
///////////////////////////

using namespace OpenMS;

START_TEST(AxisTickCalculator, "$Id: AxisTickCalculator_test.C 167 cfriedle $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CHECK((static void calcGridLines(double x1, double x2, int levels, GridVector& grid)))
	std::vector<std::vector<double> > vector;
	std::vector<double > vec1(4);
	vec1[0] = 1.0;
	vec1[1] = 2.0;
	vec1[2] = 3.0;
	vec1[3] = 4.0;
	vector.push_back(vec1);
	std::vector<double > vec2(3);
	vec2[0] = 1.5;
	vec2[1] = 2.5;
	vec2[2] = 3.5;
	vector.push_back(vec2);
	std::vector<double > vec3(6);
	vec3[0]= 1.25;
	vec3[1]= 1.75;
	vec3[2]= 2.25;
	vec3[3]= 2.75;
	vec3[4]= 3.25;
	vec3[5]= 3.75;
	vector.push_back(vec3);
	std::vector<std::vector<double> > vector1;
	AxisTickCalculator::calcGridLines(1.0,4.0,3,vector1);
	
	TEST_REAL_EQUAL(vector.size(),vector1.size());
	TEST_REAL_EQUAL(vector[0].size(),vector1[0].size());
	TEST_REAL_EQUAL(vector[1].size(),vector1[1].size());
	TEST_REAL_EQUAL(vector[2].size(),vector1[2].size());
	
	TEST_REAL_EQUAL(vector[0][0],vector1[0][0]);
	TEST_REAL_EQUAL(vector[0][1],vector1[0][1]);
	TEST_REAL_EQUAL(vector[0][2],vector1[0][2]);
	TEST_REAL_EQUAL(vector[0][3],vector1[0][3]);
	
	TEST_REAL_EQUAL(vector[1][0],vector1[1][0]);
	TEST_REAL_EQUAL(vector[1][1],vector1[1][1]);
	TEST_REAL_EQUAL(vector[1][2],vector1[1][2]);
	
	TEST_REAL_EQUAL(vector[2][0],vector1[2][0]);
	TEST_REAL_EQUAL(vector[2][1],vector1[2][1]);
	TEST_REAL_EQUAL(vector[2][2],vector1[2][2]);
	TEST_REAL_EQUAL(vector[2][3],vector1[2][3]);
	TEST_REAL_EQUAL(vector[2][4],vector1[2][4]);
	TEST_REAL_EQUAL(vector[2][5],vector1[2][5]);
RESULT


CHECK((static void calcLogGridLines(double x1, double x2, GridVector& grid)))
	std::vector<std::vector<double> > vector;
	std::vector<double > vec1(1);
	vec1[0] = 1.0;
	vector.push_back(vec1);
	std::vector<double > vec2(8);
	vec2[0] = 1.30103;
	vec2[1] = 1.47712;
	vec2[2] = 1.60206;
	vec2[3] = 1.69897;
	vec2[4] = 1.77815;
	vec2[5] = 1.8451;
	vec2[6] = 1.90309;
	vec2[7] = 1.95425;
	
	vector.push_back(vec2);
	
	std::vector<std::vector<double> > vector1;
	AxisTickCalculator::calcLogGridLines(log10(10.0),log10(100.0),vector1);
	
	TEST_REAL_EQUAL(vector.size(),vector1.size());
	
	TEST_REAL_EQUAL(vector[0].size(),vector1[0].size());
	TEST_REAL_EQUAL(vector[1].size(),vector1[1].size());
	
	TEST_REAL_EQUAL(vector[0][0],vector1[0][0]);
	
	
	TEST_REAL_EQUAL(vector[1][0],vector1[1][0]);
	TEST_REAL_EQUAL(vector[1][1],vector1[1][1]);
	TEST_REAL_EQUAL(vector[1][2],vector1[1][2]);
	TEST_REAL_EQUAL(vector[1][3],vector1[1][3]);
	TEST_REAL_EQUAL(vector[1][4],vector1[1][4]);
	TEST_REAL_EQUAL(vector[1][5],vector1[1][5]);
	TEST_REAL_EQUAL(vector[1][6],vector1[1][6]);
	TEST_REAL_EQUAL(vector[1][7],vector1[1][7]);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



