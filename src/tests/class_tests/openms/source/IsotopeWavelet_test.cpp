// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletConstants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>

using namespace OpenMS;
using namespace std;


START_TEST(IsotopeWavelet, "$Id$")

IsotopeWavelet* nullPointer = nullptr;
START_SECTION((static IsotopeWavelet* getInstance()))
  TEST_EQUAL(IsotopeWavelet::getInstance(), nullPointer)
END_SECTION

START_SECTION(static UInt getMaxCharge())
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 1)
END_SECTION

START_SECTION(static Size getGammaTableMaxIndex())
	TEST_EQUAL(IsotopeWavelet::getGammaTableMaxIndex(), 0)
END_SECTION

START_SECTION(static Size getExpTableMaxIndex())
	TEST_EQUAL(IsotopeWavelet::getExpTableMaxIndex(), 0)
END_SECTION

START_SECTION((static void setMaxCharge(const UInt max_charge))) 
	IsotopeWavelet::setMaxCharge(3);
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 3)
END_SECTION

START_SECTION((static double getTableSteps()))
	TEST_NOT_EQUAL(IsotopeWavelet::getTableSteps(), 0)
END_SECTION

START_SECTION((static void setTableSteps(const double table_steps))) 
	IsotopeWavelet::setTableSteps(0.0001);
	TEST_EQUAL(IsotopeWavelet::getTableSteps(), 0.0001)
END_SECTION

START_SECTION((static double getInvTableSteps())) 
	IsotopeWavelet::getInvTableSteps();
	TEST_EQUAL(IsotopeWavelet::getInvTableSteps(), 10000)
END_SECTION

START_SECTION((static double getLambdaL(const double m)))
	TEST_REAL_SIMILAR(IsotopeWavelet::getLambdaL(1000), 0.75632)
END_SECTION


START_SECTION((static UInt getMzPeakCutOffAtMonoPos(const double mass, const UInt z)))
	TEST_EQUAL(IsotopeWavelet::getMzPeakCutOffAtMonoPos(1000, 1), 5)
END_SECTION

START_SECTION((static UInt getNumPeakCutOff(const double mass, const UInt z)))
	TEST_EQUAL(IsotopeWavelet::getNumPeakCutOff(1000, 1), 4)
END_SECTION

START_SECTION((static UInt getNumPeakCutOff(const double mz)))
	TEST_EQUAL(IsotopeWavelet::getNumPeakCutOff(1000), 4)
END_SECTION


IsotopeWavelet* iw = nullptr;
START_SECTION((static IsotopeWavelet* init(const double max_m, const UInt max_charge)))
	iw = IsotopeWavelet::init (4000, 4);
  TEST_NOT_EQUAL(iw, nullPointer)
	TEST_EQUAL (IsotopeWavelet::getMaxCharge(), 4)
END_SECTION


UInt size=0;
START_SECTION((static const IsotopeDistribution::ContainerType& getAveragine (const double m, UInt* size=NULL)))
	IsotopeWavelet::getAveragine (1000, &size);
	TEST_EQUAL (size, 4)	 
END_SECTION


double v=-1;
START_SECTION((static double getValueByMass (const double t, const double m, const UInt z, const Int mode=+1))) 
	TOLERANCE_ABSOLUTE (1e-4)
	for (UInt c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByMass (Constants::IW_HALF_NEUTRON_MASS/(c+1.), 1000, c+1, 1);
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION

START_SECTION((static double getValueByLambda (const double lambda, const double tz1))) 
	for (Size c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByLambda (iw->getLambdaL(1000*(c+1)-(c+1)*Constants::IW_PROTON_MASS), Constants::IW_HALF_NEUTRON_MASS*(c+1)+1);
		TOLERANCE_ABSOLUTE (1e-4)
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION

START_SECTION((static double getValueByLambdaExtrapol (const double lambda, const double tz1))) 
	for (Size c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByLambdaExtrapol (iw->getLambdaL(1000*(c+1)-(c+1)*Constants::IW_PROTON_MASS), Constants::IW_HALF_NEUTRON_MASS*(c+1)+1);
		TOLERANCE_ABSOLUTE (1e-4)
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION

START_SECTION((static double getValueByLambdaExact (const double lambda, const double tz1))) 
	for (Size c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByLambdaExact (iw->getLambdaL(1000*(c+1)-(c+1)*Constants::IW_PROTON_MASS), Constants::IW_HALF_NEUTRON_MASS*(c+1)+1);
		TOLERANCE_ABSOLUTE (1e-4)
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION


START_SECTION(static float myPow(float a, float b))
	TEST_EQUAL (int(IsotopeWavelet::myPow(1.1f, 3.0f)*10), 13);
END_SECTION

START_SECTION(static void destroy())
	IsotopeWavelet::destroy();
	TEST_EQUAL (IsotopeWavelet::getExpTableMaxIndex(), 0);
END_SECTION

END_TEST
