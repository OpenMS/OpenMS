// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Witold Wolski, Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#ifdef USE_BOOST_UNIT_TEST
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>

// make work with OpenMS
#define START_TEST(val)
#define END_TEST
#define END_SECTION

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5)
#define TEST_REAL_SIMILAR(val1, val2) \
  BOOST_CHECK ( boost::test_tools::check_is_close(val1, val2, boost::test_tools::fraction_tolerance(1.e-5) ));
#include <iomanip>
#else
#include <OpenMS/CONCEPT/ClassTest.h>
#endif

using namespace std;
using namespace OpenMS;
//using namespace OpenMS::OpenSWATH;

///////////////////////////

START_TEST(DIAHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(getBYSeries_test)
{
	String sequence = "SYVAWDR";
	std::vector<double> bseries, yseries;
	OpenMS::AASequence a = OpenMS::AASequence(sequence);
  OpenMS::DIAHelpers::getBYSeries(a, bseries, yseries);
	bseries.clear();
  OpenMS::DIAHelpers::getTheorMasses(a, bseries);

}
END_SECTION

START_SECTION(getAveragineIsotopeDistribution_test)
{

	std::vector<std::pair<double, double> > tmp;
	OpenMS::DIAHelpers::getAveragineIsotopeDistribution(100., tmp);
	TEST_EQUAL(tmp.size() == 4, true);

	double mass1[] = { 100, 101.00048, 102.00096, 103.00144 };
	double int1[] =
			{ 0.9512718332, 0.04579662689, 0.002828078664, 0.0001016459634 };

	double * mm = &mass1[0];
	double * ii = &int1[0];
	for (int i = 0; i < tmp.size(); ++i, ++mm, ++ii) {

		std::cout << "mass :" << std::setprecision(10) << tmp[i].first
				<< "intensity :" << tmp[i].second << std::endl;
		TEST_REAL_SIMILAR(tmp[i].first, *mm);
		TEST_REAL_SIMILAR(tmp[i].second, *ii);
	}

	tmp.clear();
	OpenMS::DIAHelpers::getAveragineIsotopeDistribution(30., tmp);
	double mass2[] = { 30, 31.0005, 32.001, 33.0014 };
	double int2[] = { 0.989072, 0.010925, 2.4738e-06, 1.41508e-10 };
	mm = &mass2[0];
	ii = &int2[0];
	for (int i = 0; i < tmp.size(); ++i, ++mm, ++ii) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
		std::cout << "mass :" << std::setprecision(10) << tmp[i].first
				<< "intensity :" << tmp[i].second << std::endl;
		std::cout << i << "dm" <<  *mm - tmp[i].first << " di " << *ii - tmp[i].second << std::endl;
		TEST_REAL_SIMILAR(tmp[i].first, *mm)
		TEST_REAL_SIMILAR(tmp[i].second, *ii)
	}

	tmp.clear();
	OpenMS::DIAHelpers::getAveragineIsotopeDistribution(110., tmp);
	for (int i = 0; i < tmp.size(); ++i) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
	}

	tmp.clear();
	OpenMS::DIAHelpers::getAveragineIsotopeDistribution(120., tmp);
	for (int i = 0; i < tmp.size(); ++i) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
	}

	tmp.clear();
	OpenMS::DIAHelpers::getAveragineIsotopeDistribution(300., tmp);
	for (int i = 0; i < tmp.size(); ++i) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
	}

	tmp.clear();
	OpenMS::DIAHelpers::getAveragineIsotopeDistribution(500., tmp);
	for (int i = 0; i < tmp.size(); ++i) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
	}

}
END_SECTION

START_SECTION(simulateSpectrumFromAASequence_test)
{
	String sequence = "SYVAWDR";
	OpenMS::AASequence a = OpenMS::AASequence(sequence);
	std::vector<double> masses1;
	std::vector<std::pair<double, double> > tmp, out;
	OpenMS::DIAHelpers::simulateSpectrumFromAASequence(a, masses1, tmp);

	std::copy(masses1.begin(), masses1.end(),
			std::ostream_iterator<double>(std::cout, " "));
	std::cout << std::endl;
	for (int i = 0; i < tmp.size(); ++i) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
	}
	OpenMS::DIAHelpers::modifyMassesByCharge(tmp, out, 2.);
	OpenMS::DIAHelpers::addPreisotopeWeights(masses1, tmp);
	std::cout << "preisotope weights added" << std::endl;

	for (int i = 0; i < tmp.size(); ++i) {
		std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
				<< std::endl;
	}

}
END_SECTION

START_SECTION(addIsotopesToSpec_test)
{
	std::vector<std::pair<double, double> > tmp_, out;
	tmp_.push_back(std::make_pair(100., 100.));
	tmp_.push_back(std::make_pair(200., 300.));
	tmp_.push_back(std::make_pair(300., 200.));

  OpenMS::DIAHelpers::addIsotopes2Spec(tmp_, out);
	std::cout << "addIsotopesToSpec_test" << std::endl;
	for (int i = 0; i < out.size(); ++i) {
		std::cout << out[i].first << " " << out[i].second << std::endl;
	}

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
