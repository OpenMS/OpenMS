// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/MassExplainer.h>

using namespace OpenMS;
using namespace std;

START_TEST(ILPDCWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ILPDCWrapper* ptr = 0;
ILPDCWrapper* null_ptr = 0;
START_SECTION(ILPDCWrapper())
{
	ptr = new ILPDCWrapper();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~ILPDCWrapper())
{
	delete ptr;
}
END_SECTION


START_SECTION((DoubleReal compute(const FeatureMap<> fm, PairsType &pairs, Size verbose_level) const))
{
  EmpiricalFormula ef("H1");
  Adduct a(+1, 1, ef.getMonoWeight(), "H1", 0.1, 0, "");
  MassExplainer::AdductsType potential_adducts_;
  potential_adducts_.push_back(a);
  MassExplainer me(potential_adducts_, 1, 3, 2, 0, 0);
  FeatureMap<> fm;
  ILPDCWrapper::PairsType pairs;

  ILPDCWrapper iw;
  iw.compute(fm, pairs, 1);

  // check that it runs without pairs (i.e. all clusters are singletons)
  TEST_EQUAL(pairs.size(), 0);

  // real data test
  

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



