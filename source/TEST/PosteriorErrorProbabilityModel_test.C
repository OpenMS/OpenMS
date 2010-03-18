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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <vector>
#include <iostream>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(PosteriorErrorProbabilityModel, "$Id:$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PosteriorErrorProbabilityModel* ptr = 0;
START_SECTION(PosteriorErrorProbabilityModel())
{
	ptr = new PosteriorErrorProbabilityModel();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((virtual ~PosteriorErrorProbabilityModel()))
{
  delete ptr;
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void fit(std::vector<double>& x_scores, std::vector<double>& probabilites)))
	ptr = new PosteriorErrorProbabilityModel();
	vector<double> score_vector;
	
	score_vector.push_back(-0.39);
	score_vector.push_back(0.06);
	score_vector.push_back(0.12);
	score_vector.push_back(0.48);
	score_vector.push_back(0.94);
	score_vector.push_back(1.01);
	score_vector.push_back(1.67);
	score_vector.push_back(1.68);
	score_vector.push_back(1.76);
	score_vector.push_back(1.80);
	score_vector.push_back(2.44);
	score_vector.push_back(3.25);
	score_vector.push_back(3.72);
	score_vector.push_back(4.12);
	score_vector.push_back(4.28);
	score_vector.push_back(4.60);
	score_vector.push_back(4.92);
	score_vector.push_back(5.28);
	score_vector.push_back(5.53);
	score_vector.push_back(6.22);
	
	
	vector<double> probabilities;
	Param param;
	param.setValue("number_of_bins", 10);
	ptr->setParameters(param);
	ptr->fit(score_vector, probabilities);
	
	Size i(0),j(1);
	TOLERANCE_ABSOLUTE(0.5)
	TEST_REAL_SIMILAR(ptr->getGaussFitResult().x0 , 4.62)
	TEST_REAL_SIMILAR(ptr->getGaussFitResult().sigma, 0.87)
	TEST_REAL_SIMILAR(ptr->getGumbelFitResult().a, 1.06)
	TEST_REAL_SIMILAR(ptr->getGumbelFitResult().b, 0.77)
	TEST_REAL_SIMILAR(ptr->getNegativePrior(), 0.546)
	TOLERANCE_ABSOLUTE(0.001)
	while(i < score_vector.size() && j < score_vector.size())
	{
		cout<<"i: "<<score_vector[i] << ", j: "<<score_vector[j]<<endl;
		cout<<"pi:"<<probabilities[i] <<", j: "<<probabilities[j]<<endl;
		if(score_vector[i] <= score_vector[j])
		{
			TEST_EQUAL(probabilities[i] >= probabilities[j],true)
		}
		else
		{
			TEST_EQUAL(probabilities[i] >= probabilities[j],true)
		}
	++i;
	++j;
	}
END_SECTION

START_SECTION((GaussFitter::GaussFitResult getGaussFitResult() const))
//tested in fit
NOT_TESTABLE
END_SECTION

START_SECTION((GumbelDistributionFitter::GumbelDistributionFitResult getGumbelFitResult() const))
//tested in fit
NOT_TESTABLE
END_SECTION

START_SECTION((DoubleReal getNegativePrior() const))
//tested in fit
NOT_TESTABLE
END_SECTION

START_SECTION((const String getGumbelGnuplotFormula() const))
String gumbel = ptr->getGumbelGnuplotFormula();
//"f(x)= = (1/0.907832") * exp(( 1.48185 - x)/0.907832) * exp(-exp(( 1.48185 - x)/0.907832))"

	TEST_EQUAL(gumbel.hasSubstring("f(x)="), true)
	TEST_EQUAL(gumbel.hasSubstring("(1/0.907832)"), true)
	TEST_EQUAL(gumbel.hasSubstring("exp(( 1.48185- x)/0.907832)"), true)
	TEST_EQUAL(gumbel.hasSubstring(") * exp(-exp(("), true)
END_SECTION
				
START_SECTION((const String getGaussGnuplotFormula() const))
String gauss = ptr->getGaussGnuplotFormula();
//g(x)=0.444131 * exp(-(x - 5.05539) ** 2 / 2 / (0.898253) ** 2)
	TEST_EQUAL(gauss.hasSubstring("g(x)="), true)
	TEST_EQUAL(gauss.hasSubstring(" * exp(-(x - "), true)
	TEST_EQUAL(gauss.hasSubstring(") ** 2 / 2 / ("), true)
	TEST_EQUAL(gauss.hasSubstring(") ** 2)"), true)
END_SECTION

START_SECTION((const String getBothGnuplotFormula(double negative_prior) const))
NOT_TESTABLE
delete ptr;
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



