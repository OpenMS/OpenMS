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
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DIAHelpers.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/timer.hpp>

#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5)
#define END_SECTION
#define TEST_REAL_SIMILAR(val1, val2) \
  BOOST_CHECK ( boost::test_tools::check_is_close(val1, val2, boost::test_tools::fraction_tolerance(1.e-5) ));
#define TEST_EQUAL(val1, val2) BOOST_CHECK_EQUAL(val1, val2);

using namespace std;

BOOST_AUTO_TEST_CASE(testIntegrateWindows_test)
{


	OpenSwath::SpectrumPtr spec(new OpenSwath::Spectrum());
	OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);

	mass->data.push_back(100.);
	mass->data.push_back(101.);
	mass->data.push_back(102.);
	mass->data.push_back(103.);
	mass->data.push_back(104.);
	mass->data.push_back(105.);
	mass->data.push_back(106.);

	OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
	intensity->data.push_back(1.);
	intensity->data.push_back(1.);
	intensity->data.push_back(1.);
	intensity->data.push_back(1.);
	intensity->data.push_back(1.);
	intensity->data.push_back(1.);
	intensity->data.push_back(1.);
	spec->binaryDataArrayPtrs.push_back(mass);
	spec->binaryDataArrayPtrs.push_back(intensity);
	double mz, intens;
	OpenSwath::integrateWindow(spec, 101., 103., mz, intens);
	std::cout << "mz : " << mz << " int : " << intens << std::endl;
	std::vector<double> windows, intInt, intMz;
	windows.push_back(101.);
	windows.push_back(103.);
	windows.push_back(105.);
	OpenSwath: integrateWindows(spec, windows, 2, intInt, intMz);

	std::cout << "print Int" << std::endl;
	std::copy(intInt.begin(), intInt.end(),
			std::ostream_iterator<double>(std::cout, " "));
	std::cout << std::endl << "print mz" << intMz.size() << std::endl;
	std::cout << intMz[0] << " " << intMz[1] << " " << intMz[2] << std::endl;
	std::copy(intMz.begin(), intMz.end(),
			std::ostream_iterator<double>(std::cout, " "));

}
END_SECTION

BOOST_AUTO_TEST_CASE(testDotProdScore)
{
	double arr1[] = { 100., 200., 4., 30., 20. };
	double arr2[] = { 100., 100., 4., 100., 200. };
	std::vector<double> vec1;
	std::vector<double> vec2;
	vec1.assign(arr1, arr1 + sizeof(arr1) / sizeof(double));
	vec2.assign(arr2, arr2 + sizeof(arr2) / sizeof(double));
	/*
	x<-c(100., 200., 4., 30., 20.)
	y<-c(100., 100., 4., 100., 200.)
	xs<-sqrt(x)
	ys<-sqrt(y)
	xsn<-xs/sqrt(sum(xs*xs))
	ysn<-ys/sqrt(sum(ys*ys))
	sum(xsn*ysn)
	*/
	//0.8604286
	double scor = OpenSwath::dotprodScoring(vec1,vec2);

	BOOST_CHECK ( boost::test_tools::check_is_close(scor, 0.8604286, boost::test_tools::fraction_tolerance(1.e-5) ));

	//xsm <- xs/sum(xs)
	//ysm <-ys/sum(ys)
	//sum(fabs(ysm-xsm))
	scor = OpenSwath::manhattanScoring(vec1,vec2);
	BOOST_CHECK ( boost::test_tools::check_is_close(scor, 0.4950837, boost::test_tools::fraction_tolerance(1.e-5) ));
	//0.4950837
}
END_SECTION
