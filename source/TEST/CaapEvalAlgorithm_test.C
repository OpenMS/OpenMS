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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithmRecall.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class CEA
	 : public CaapEvalAlgorithm
	{
		public:
			void evaluate(const ConsensusMap&, const ConsensusMap&, DoubleReal& real)
			{
				real = 1.5;
			}
	};
}

START_TEST(CaapEval, "$Id CaapEvalAlgorithm_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CEA* ptr = 0;
START_SECTION((CaapEvalAlgorithm()))
	ptr = new CEA();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~CaapEvalAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap& mapin1, const ConsensusMap& mapin2, DoubleReal& realin)=0))
	CEA cea;
	ConsensusMap map1;
	ConsensusMap map2;
	DoubleReal real;
	cea.evaluate(map1, map2, real);
	TEST_EQUAL(real, 1.5)
END_SECTION

//isSameHandle testen!?

START_SECTION((static void registerChildren()))
{
	TEST_STRING_EQUAL(Factory<CaapEvalAlgorithm>::registeredProducts()[0],CaapEvalAlgorithmPrecision::getProductName());
	TEST_STRING_EQUAL(Factory<CaapEvalAlgorithm>::registeredProducts()[1],CaapEvalAlgorithmRecall::getProductName());
	TEST_EQUAL(Factory<CaapEvalAlgorithm>::registeredProducts().size(),2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
